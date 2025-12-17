import uproot
import numpy as np
import matplotlib.pyplot as plt
import os
import awkward as ak

# CONSTANTS
SOURCE_ENERGIES_KEV = [200, 2000, 10000]
ORIENTATIONS = [0, 90, 180]

TREE_KEYS = ["SToGS;1", "SToGS;2"]

BIN_MAP = { # used for histogram displays
    200:  np.linspace(0, 0.2, 500),
    2000: np.linspace(0, 2, 500),
    10000: np.linspace(0, 11, 500),
}

fig, ax = plt.subplots(len(SOURCE_ENERGIES_KEV), len(ORIENTATIONS))
fig.subplots_adjust(hspace=0.4)


# DATA
for orient_idx, orientation in enumerate(ORIENTATIONS):
    for row_idx, energy_kev in enumerate(SOURCE_ENERGIES_KEV): # one plot (two detectors) per iteration

        # FILE NAMING DEFINITION
        filename = f"two_tubes//two_{orientation}_{energy_kev}.root"
        print(f"Processing file: {filename}...")

        if not os.path.exists(filename):
            print(f"WARNING: File not found: {filename}. Skipping this energy.")
            continue
        
        for tree_key in TREE_KEYS:

            # arraye do energi na jeden symulacje
            sum_EvE_0 = []
            sum_EvE_1 = []

            print(f"Processing tree: {tree_key}...")

            tree = uproot.open(filename)[tree_key]
            EvE   = tree["Ev.E"].array()     # jagged array: shape (events, hits)
            EvUID = tree["Ev.UID"].array()   # jagged array: same structure

            #------experoent----- # works
            has0 = ak.any(EvUID == 0, axis=1)
            has1 = ak.any(EvUID == 1, axis=1)

            coinc_mask = has0 & has1

            EvE   = EvE[coinc_mask]
            EvUID = EvUID[coinc_mask]
            #--------------------

            # Create event-level masks for UID == 0 and UID == 1
            mask0 = EvUID == 0
            mask1 = EvUID == 1

            # Apply the masks and compute per-event sums (vectorized)
            sum_EvE_0 = ak.sum(EvE[mask0], axis=1)
            sum_EvE_1 = ak.sum(EvE[mask1], axis=1)

            # deleting NON coincidence
            mask = (sum_EvE_0 != 0) & (sum_EvE_1 != 0)

            sum_EvE_0 = sum_EvE_0[mask]
            sum_EvE_1 = sum_EvE_1[mask]
            # n o coincidence
            n_co = ak.sum(mask) 

        # Plotting the histograms
        ax[row_idx,orient_idx].scatter(sum_EvE_0, sum_EvE_1, color='green', s=5)
        ax[row_idx,orient_idx].set_title(f"{energy_kev / 1000} MeV, Relative orientation: {orientation}, number of coincidence: {n_co}", fontsize=10)
        ax[row_idx,orient_idx].grid(True, alpha=0.5)
        ax[row_idx,orient_idx].set_ylabel("E1")
        ax[row_idx,orient_idx].set_xlabel("E0")       
        


fig.suptitle(f"Coincidence", fontsize=14, y=0.98)
print("Processing done, witing for plots...")
# --- 5. Final Display ---

plt.show()