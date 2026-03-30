/*
========================================================================================
    POPFINDER CLASSIFICATION (Birchard et al.)
========================================================================================
    Uses the popfinder Python package (Birchard, Boccia, Lounder, Colston-Nepali
    & Friesen) for neural network-based genetic population assignment.

    PopFinder trains a classifier neural network (PyTorch) with k-fold
    cross-validation for accuracy estimation, then assigns unknown samples.

    Input genotype TSV matrices are converted to HDF5 format for popfinder
    compatibility. Dummy coordinates (x=0, y=0) are used since the classifier
    mode does not require geographic data.

    Only runs at nodes where the majority of classes have sufficient samples.
*/

process RUN_POPFINDER {
    tag "${node_dir.name}"
    label 'process_high'

    publishDir "${params.outdir}/classification/popfinder", mode: 'copy'

    input:
    tuple path(node_dir), path(node_meta)

    output:
    tuple val("${node_dir.name}"), path("${node_dir.name}_popfinder_results"), emit: results

    script:
    """
    #!/usr/bin/env python3

    import json
    import os
    import sys
    import numpy as np
    import pandas as pd
    import h5py
    from pathlib import Path

    node_name = os.path.basename("${node_dir}")
    out_dir = f"{node_name}_popfinder_results"
    os.makedirs(out_dir, exist_ok=True)

    with open("${node_meta}") as f:
        meta = json.load(f)

    print(f"popfinder at node: {node_name}")
    print(f"Groups: {meta['groups']}")

    # Load reference data
    ref_geno = pd.read_csv(os.path.join("${node_dir}", "reference_genotypes.tsv"), sep="\\t")
    ref_labels = pd.read_csv(os.path.join("${node_dir}", "reference_labels.tsv"), sep="\\t")

    sample_ids = ref_geno["sample_id"].values
    X_ref = ref_geno.drop("sample_id", axis=1).values.astype(float)
    labels = ref_labels.set_index("sample_id").loc[sample_ids, "group"].values

    # Impute missing with column means
    col_means = np.nanmean(X_ref, axis=0)
    nan_idx = np.where(np.isnan(X_ref))
    if len(nan_idx[0]) > 0:
        X_ref[nan_idx] = np.take(col_means, nan_idx[1])

    unique_labels = sorted(set(labels))
    n_samples, n_snps = X_ref.shape
    n_classes = len(unique_labels)
    print(f"Data: {n_samples} samples x {n_snps} SNPs, {n_classes} classes")

    # --- Handle single-group case ---
    if n_classes < 2:
        single_grp = unique_labels[0]
        print(f"Only 1 group -- skipping NN, assigning all to: {single_grp}")

        loocv_df = pd.DataFrame({
            "sample_id": sample_ids, "true_group": labels,
            "predicted_group": labels, "posterior_max": 1.0, "correct": True
        })
        loocv_df.to_csv(os.path.join(out_dir, "loocv_results.tsv"), sep="\\t", index=False)

        unk_file = os.path.join("${node_dir}", "unknown_genotypes.tsv")
        if os.path.exists(unk_file):
            unk_geno = pd.read_csv(unk_file, sep="\\t")
            if len(unk_geno) > 0:
                unk_df = pd.DataFrame({
                    "sample_id": unk_geno["sample_id"].values,
                    "predicted_group": single_grp,
                    "probability": 1.0
                })
                unk_df[f"posterior_{single_grp}"] = 1.0
                unk_df.to_csv(os.path.join(out_dir, "unknown_predictions.tsv"), sep="\\t", index=False)

        summary_df = pd.DataFrame([{
            "node": node_name, "method": "popfinder",
            "n_reference": n_samples, "n_snps": n_snps, "n_groups": 1,
            "loocv_accuracy": 1.0, "class_imbalance_ratio": 1.0
        }])
        summary_df.to_csv(os.path.join(out_dir, "summary.tsv"), sep="\\t", index=False)
        print(f"popfinder complete for {node_name} (single group)")
        sys.exit(0)

    # --- Load unknown genotypes if present ---
    unk_ids = None
    X_unk = None
    unk_file = os.path.join("${node_dir}", "unknown_genotypes.tsv")
    if os.path.exists(unk_file):
        unk_geno = pd.read_csv(unk_file, sep="\\t")
        if len(unk_geno) > 0:
            unk_ids = unk_geno["sample_id"].values
            X_unk = unk_geno.drop("sample_id", axis=1).values.astype(float)
            # Impute unknowns with reference means
            unk_nan = np.where(np.isnan(X_unk))
            if len(unk_nan[0]) > 0:
                X_unk[unk_nan] = np.take(col_means, unk_nan[1])

    # --- Convert to HDF5 for popfinder ---
    # popfinder expects HDF5 with 'derived_counts' (samples x SNPs) and 'samples'
    all_ids = list(sample_ids)
    all_X = X_ref.copy()
    all_labels = list(labels)

    if unk_ids is not None:
        all_ids.extend(list(unk_ids))
        all_X = np.vstack([X_ref, X_unk])
        all_labels.extend([None] * len(unk_ids))  # unknowns have no pop

    hdf5_path = "genotypes.hdf5"
    with h5py.File(hdf5_path, "w") as hf:
        hf.create_dataset("derived_counts", data=all_X.astype(int))
        hf.create_dataset("samples", data=np.array(all_ids, dtype=h5py.string_dtype()))

    # --- Create sample_data TSV ---
    # popfinder requires: sampleID, pop, x, y
    # Classifier mode doesn't use x/y but columns must exist
    sample_data_path = "sample_data.tsv"
    sample_df = pd.DataFrame({
        "sampleID": all_ids,
        "pop": [l if l is not None else np.nan for l in all_labels],
        "x": 0.0,
        "y": 0.0
    })
    sample_df.to_csv(sample_data_path, sep="\\t", index=False)

    # --- Run popfinder ---
    import torch
    # popfinder 0.0.6 calls torch.load() without weights_only=False,
    # which fails on PyTorch 2.6+ (default changed to weights_only=True).
    # Monkey-patch torch.load to restore old behavior for popfinder's own models.
    _orig_torch_load = torch.load
    def _torch_load_compat(*args, **kwargs):
        kwargs.setdefault("weights_only", False)
        return _orig_torch_load(*args, **kwargs)
    torch.load = _torch_load_compat

    from popfinder.dataloader import GeneticData
    from popfinder.classifier import PopClassifier

    # Monkey-patch _sort_samples for pandas 2.x compatibility
    # Two fixes: (1) decode bytes from HDF5 so reindex matches string TSV index,
    # (2) use .iloc[x] instead of [x] for positional access after set_index.
    def _sort_samples_fixed(self, locs, samples):
        if not pd.Series(["x", "pop", "y", "sampleID"]).isin(locs.columns).all():
            raise ValueError("sample_data does not have correct columns")
        # Decode bytes to str if HDF5 stored as bytes
        samples = [s.decode() if isinstance(s, bytes) else str(s) for s in samples]
        locs["id"] = locs["sampleID"]
        locs.set_index("id", inplace=True)
        locs = locs.reindex(samples)
        locs["order"] = np.arange(0, len(locs))
        if not all(
            [locs["sampleID"].iloc[x] == samples[x] for x in range(len(samples))]
        ):
            raise ValueError("sample ordering failed! Check that sample IDs match VCF.")
        return locs
    GeneticData._sort_samples = _sort_samples_fixed

    print("Initializing popfinder GeneticData...")
    data = GeneticData(
        genetic_data=hdf5_path,
        sample_data=sample_data_path,
        test_size=0.2,
        seed=42
    )

    print(f"Known samples: {len(data.knowns)}, Unknown: {len(data.unknowns)}")

    # Initialize classifier
    pf_out = os.path.join(out_dir, "popfinder_output")
    classifier = PopClassifier(data=data, random_state=42, output_folder=pf_out)

    # Class imbalance check
    class_counts = data.knowns["pop"].value_counts()
    max_ratio = class_counts.max() / class_counts.min()
    print(f"Class sizes: {dict(class_counts)}")
    print(f"Max class ratio: {max_ratio:.1f}")

    # Pre-fit label encoder on ALL classes so CV folds with rare classes don't fail
    from sklearn.preprocessing import LabelEncoder
    all_pops = sorted(data.knowns["pop"].unique())
    classifier.label_enc = LabelEncoder()
    classifier.label_enc.fit(all_pops)

    # Monkey-patch _split_input_classifier: the original re-fits label_enc on
    # just the training fold via fit_transform(), which drops rare classes.
    # Fix: use transform() only (encoder already fit on all classes above).
    from popfinder._helper import _data_converter
    import popfinder._helper as _pf_helper
    import popfinder.classifier as _pf_classifier
    def _split_fixed(clf, input):
        train_input, valid_input = input
        X_train = train_input["alleles"]
        X_valid = valid_input["alleles"]
        y_train = clf.label_enc.transform(train_input["pop"])
        y_valid = clf.label_enc.transform(valid_input["pop"])
        X_train, y_train = _data_converter(X_train, y_train)
        X_valid, y_valid = _data_converter(X_valid, y_valid)
        return X_train, y_train, X_valid, y_valid
    # Patch in both modules (classifier.py has its own local reference)
    _pf_helper._split_input_classifier = _split_fixed
    _pf_classifier._split_input_classifier = _split_fixed

    # Adaptive CV splits: min class size must be sufficient for stratified k-fold + test holdout
    min_class_n = class_counts.min()
    cv_splits = min(5, int(min_class_n))
    if min_class_n < 5:
        print(f"WARNING: smallest class has {min_class_n} samples, too few for neural net CV. Skipping popfinder.")
        # Write minimal output so downstream processes handle NO_POPFINDER
        summary_df = pd.DataFrame([{"node": node_name, "loocv_accuracy": np.nan,
                                     "test_accuracy": np.nan, "note": "skipped_small_class"}])
        summary_df.to_csv(os.path.join(out_dir, "summary.tsv"), sep="\\t", index=False)
        print(f"popfinder complete for {node_name} (skipped)")
        sys.exit(0)
    print(f"Training popfinder classifier with {cv_splits}-fold CV...")
    classifier.train(
        epochs=100,
        valid_size=0.2,
        cv_splits=cv_splits,
        nreps=1,
        learning_rate=0.001,
        batch_size=32,
        dropout_prop=0.25,
        overwrite_results=True
    )

    # Get CV accuracy (popfinder 0.0.6 bug: __cv_accuracy is never set)
    try:
        cv_accuracy = classifier.cv_accuracy
        if cv_accuracy is not None:
            print(f"popfinder CV accuracy: {cv_accuracy * 100:.1f}%")
        else:
            cv_accuracy = 0.0
    except AttributeError:
        cv_accuracy = 0.0
        print("Warning: CV accuracy attribute not set (popfinder bug)")

    # Test on held-out set
    try:
        classifier.test(use_best_model=True)
        test_results = classifier.test_results
        test_accuracy = classifier.accuracy
        print(f"popfinder test accuracy: {test_accuracy * 100:.1f}%")
    except AttributeError:
        test_accuracy = 0.0
        print("Warning: test accuracy attribute not set")

    # --- Also run LOOCV for comparability with other classifiers ---
    # popfinder's built-in CV is k-fold, not LOO. Run manual LOOCV for
    # direct comparison with DAPC and assignPOP LOOCV results.
    print("Running manual LOOCV for comparability...")

    from sklearn.preprocessing import LabelEncoder, StandardScaler
    import torch
    import torch.nn as nn
    from torch.utils.data import DataLoader, TensorDataset

    le = LabelEncoder()
    y_encoded = le.fit_transform(labels)

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_ref)

    loocv_preds = np.zeros(n_samples, dtype=int)
    loocv_probs = np.zeros((n_samples, n_classes))

    for i in range(n_samples):
        # Leave one out
        X_train = np.delete(X_scaled, i, axis=0)
        y_train = np.delete(y_encoded, i)
        X_test = X_scaled[i:i+1]

        # Simple 2-layer network matching popfinder architecture
        input_dim = X_train.shape[1]
        model = nn.Sequential(
            nn.Linear(input_dim, 128),
            nn.ReLU(),
            nn.Dropout(0.25),
            nn.Linear(128, 64),
            nn.ReLU(),
            nn.Dropout(0.25),
            nn.Linear(64, n_classes)
        )

        criterion = nn.CrossEntropyLoss()
        optimizer = torch.optim.Adam(model.parameters(), lr=0.001)

        train_X = torch.FloatTensor(X_train)
        train_y = torch.LongTensor(y_train)
        dataset = TensorDataset(train_X, train_y)
        loader = DataLoader(dataset, batch_size=32, shuffle=True)

        model.train()
        for epoch in range(100):
            for batch_X, batch_y in loader:
                optimizer.zero_grad()
                output = model(batch_X)
                loss = criterion(output, batch_y)
                loss.backward()
                optimizer.step()

        model.eval()
        with torch.no_grad():
            test_output = model(torch.FloatTensor(X_test))
            probs = torch.softmax(test_output, dim=1).numpy()[0]
            pred = np.argmax(probs)

        loocv_preds[i] = pred
        loocv_probs[i] = probs

    loocv_accuracy = np.mean(loocv_preds == y_encoded)
    print(f"LOOCV accuracy: {loocv_accuracy * 100:.1f}%")

    # Save LOOCV results
    loocv_df = pd.DataFrame({
        "sample_id": sample_ids,
        "true_group": labels,
        "predicted_group": le.inverse_transform(loocv_preds),
        "posterior_max": loocv_probs.max(axis=1),
        "correct": loocv_preds == y_encoded
    })
    loocv_df.to_csv(os.path.join(out_dir, "loocv_results.tsv"), sep="\\t", index=False)

    # --- Assign unknowns ---
    if unk_ids is not None and len(unk_ids) > 0:
        classification = None
        try:
            classifier.assign_unknown(use_best_model=True)
            # assign_unknown stores results in data.unknowns["assigned_pop"]
            unk_data = classifier.data.unknowns
            if unk_data is not None and "assigned_pop" in unk_data.columns:
                unk_results = unk_data[unk_data["sampleID"].isin(unk_ids)]
                unk_df = pd.DataFrame({
                    "sample_id": unk_results["sampleID"].values,
                    "predicted_group": unk_results["assigned_pop"].values,
                    "probability": 1.0
                })
                for label in unique_labels:
                    unk_df[f"posterior_{label}"] = np.nan
                unk_df.to_csv(os.path.join(out_dir, "unknown_predictions.tsv"), sep="\\t", index=False)
                classification = unk_df  # flag that we succeeded
        except Exception as e:
            print(f"popfinder assign_unknown failed: {e}, using LOOCV model fallback...")
            classification = None

        if classification is None:
            # Fallback: use LOOCV model for unknowns
            print("popfinder assign_unknown returned no results, using LOOCV model...")
            X_unk_scaled = scaler.transform(X_unk)

            # Train final model on all data
            final_model = nn.Sequential(
                nn.Linear(n_snps, 128), nn.ReLU(), nn.Dropout(0.25),
                nn.Linear(128, 64), nn.ReLU(), nn.Dropout(0.25),
                nn.Linear(64, n_classes)
            )
            criterion = nn.CrossEntropyLoss()
            optimizer = torch.optim.Adam(final_model.parameters(), lr=0.001)

            train_X = torch.FloatTensor(X_scaled)
            train_y = torch.LongTensor(y_encoded)
            dataset = TensorDataset(train_X, train_y)
            loader = DataLoader(dataset, batch_size=32, shuffle=True)

            final_model.train()
            for epoch in range(100):
                for batch_X, batch_y in loader:
                    optimizer.zero_grad()
                    output = final_model(batch_X)
                    loss = criterion(output, batch_y)
                    loss.backward()
                    optimizer.step()

            final_model.eval()
            with torch.no_grad():
                unk_output = final_model(torch.FloatTensor(X_unk_scaled))
                unk_probs = torch.softmax(unk_output, dim=1).numpy()
                unk_pred_idx = np.argmax(unk_probs, axis=1)

            unk_df = pd.DataFrame({
                "sample_id": unk_ids,
                "predicted_group": le.inverse_transform(unk_pred_idx),
                "probability": unk_probs.max(axis=1)
            })
            for j, label in enumerate(le.classes_):
                unk_df[f"posterior_{label}"] = unk_probs[:, j]

            unk_df.to_csv(os.path.join(out_dir, "unknown_predictions.tsv"), sep="\\t", index=False)

    # --- Save NN weights for Shiny app (pure R forward pass, no Python needed) ---
    import json as json_mod

    print("Training final model on all data for Shiny export...")
    shiny_model = nn.Sequential(
        nn.Linear(n_snps, 128), nn.ReLU(), nn.Dropout(0.25),
        nn.Linear(128, 64), nn.ReLU(), nn.Dropout(0.25),
        nn.Linear(64, n_classes)
    )
    shiny_criterion = nn.CrossEntropyLoss()
    shiny_optimizer = torch.optim.Adam(shiny_model.parameters(), lr=0.001)

    shiny_train_X = torch.FloatTensor(X_scaled)
    shiny_train_y = torch.LongTensor(y_encoded)
    shiny_dataset = TensorDataset(shiny_train_X, shiny_train_y)
    shiny_loader = DataLoader(shiny_dataset, batch_size=32, shuffle=True)

    torch.manual_seed(42)
    shiny_model.train()
    for epoch in range(100):
        for batch_X, batch_y in shiny_loader:
            shiny_optimizer.zero_grad()
            output = shiny_model(batch_X)
            loss = shiny_criterion(output, batch_y)
            loss.backward()
            shiny_optimizer.step()

    shiny_model.eval()
    state = shiny_model.state_dict()
    nn_weights = {
        "W1": state["0.weight"].numpy().tolist(),
        "b1": state["0.bias"].numpy().tolist(),
        "W2": state["3.weight"].numpy().tolist(),
        "b2": state["3.bias"].numpy().tolist(),
        "W3": state["6.weight"].numpy().tolist(),
        "b3": state["6.bias"].numpy().tolist(),
        "scaler_mean": scaler.mean_.tolist(),
        "scaler_scale": scaler.scale_.tolist(),
        "classes": le.classes_.tolist(),
        "snp_names": list(ref_geno.drop("sample_id", axis=1).columns)
    }
    nn_weights_path = os.path.join(out_dir, "nn_weights.json")
    with open(nn_weights_path, "w") as f:
        json_mod.dump(nn_weights, f)
    print(f"Saved NN weights to {nn_weights_path}")

    # --- Summary ---
    summary_df = pd.DataFrame([{
        "node": node_name,
        "method": "popfinder",
        "n_reference": n_samples,
        "n_snps": n_snps,
        "n_groups": n_classes,
        "loocv_accuracy": round(loocv_accuracy, 4),
        "kfold_cv_accuracy": round(float(cv_accuracy), 4) if cv_accuracy else 0.0,
        "test_accuracy": round(float(test_accuracy), 4) if test_accuracy else 0.0,
        "class_imbalance_ratio": round(max_ratio, 2)
    }])
    summary_df.to_csv(os.path.join(out_dir, "summary.tsv"), sep="\\t", index=False)

    print(f"popfinder complete for {node_name}")
    """
}
