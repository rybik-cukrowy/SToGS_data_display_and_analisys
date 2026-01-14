import uproot
import numpy as np
import matplotlib.pyplot as plt
import awkward as ak

# CONSTANTS
SOURCE_ENERGIES_KEV = ["multiple_gammas", "10000"]
TREE_KEYS = ["SToGS;1", "SToGS;2"]
N_DETS = 56

fig, ax = plt.subplots(1, len(SOURCE_ENERGIES_KEV), figsize=(15, 5), squeeze=False)
ax = ax[0]

for row_idx, energy_kev in enumerate(SOURCE_ENERGIES_KEV):
    coinc_map = np.zeros((N_DETS, N_DETS), dtype=int)
    filename = f"paris1/paris1_{energy_kev}.root"
    print(f"Processing file {filename}...")

    with uproot.open(filename) as f:
        for tree_key in TREE_KEYS:
            EvUID = f[tree_key]["Ev.UID"].array()

            # Flatten jagged array
            counts = ak.num(EvUID)
            event_idx_flat = np.repeat(np.arange(len(EvUID)), counts)
            detector_ids_flat = ak.to_numpy(ak.flatten(EvUID))

            # Remove duplicates per event
            combined = np.stack([event_idx_flat, detector_ids_flat], axis=1)
            unique_combined = np.unique(combined, axis=0)

            # --- Fully vectorized coincidence accumulation ---
            # First, group detector IDs by event
            # Find event boundaries
            sorted_idx = np.argsort(unique_combined[:,0])
            unique_combined = unique_combined[sorted_idx]
            event_ids = unique_combined[:,0]
            detector_ids = unique_combined[:,1]

            # Compute counts per event
            _, counts_per_event = np.unique(event_ids, return_counts=True)
            # Indices for splitting flattened array into events
            split_idx = np.cumsum(counts_per_event)

            # Split detector_ids into ragged array of events
            events_detectors = np.split(detector_ids, split_idx[:-1])

            # Vectorize all pairs across all events
            all_i = np.concatenate([np.repeat(d, len(d)) for d in events_detectors if len(d) > 1])
            all_j = np.concatenate([np.tile(d, len(d)) for d in events_detectors if len(d) > 1])

            # Accumulate into coincidence map
            np.add.at(coinc_map, (all_i, all_j), 1)

    # Remove self-coincidences
    np.fill_diagonal(coinc_map, 0)

    # Plot coincidence map
    im = ax[row_idx].imshow(coinc_map, cmap="turbo", origin="lower")
    ax[row_idx].set_title(f"{energy_kev} keV")
    ax[row_idx].set_xlabel("Detector ID")
    ax[row_idx].set_ylabel("Detector ID")
    fig.colorbar(im, ax=ax[row_idx])

print("Processing done, waiting for plots...")
plt.show()
