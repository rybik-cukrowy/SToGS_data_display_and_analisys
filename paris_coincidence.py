import uproot
import numpy as np
import matplotlib.pyplot as plt
import os
import awkward as ak

# does the plots compering coincidence between parirs of detectors

# CONSTANTS
SOURCE_ENERGIES_KEV = ["multiple_gammas"]
ID_MAIN = 0 # main detector for comarison
IDS = [1,28] # list of ids to compare to

TREE_KEYS = ["SToGS;1", "SToGS;2"]

fig, ax = plt.subplots(len(SOURCE_ENERGIES_KEV), len(IDS))
fig.subplots_adjust(hspace=0.4)

# DATA
for det_idx, det_id in enumerate(IDS):
    for row_idx, energy_kev in enumerate(SOURCE_ENERGIES_KEV): # one plot (two detectors) per iteration

        # FILE NAMING DEFINITION
        filename = f"paris1//paris1_{energy_kev}.root"
        print(f"Processing file: {filename}...")

        if not os.path.exists(filename):
            print(f"WARNING: File not found: {filename}. Skipping this energy.")
            continue
        
        for tree_key in TREE_KEYS:

            # arraye do energi na jeden porownanie i symulacje
            sum_EvE_0 = []
            sum_EvE_1 = []

            print(f"Processing tree: {tree_key}...")

            tree = uproot.open(filename)[tree_key]
            EvE   = tree["Ev.E"].array()     # jagged array: shape (events, hits)
            EvUID = tree["Ev.UID"].array()   # jagged array: same structure

            #------experoent----- # works
            has0 = ak.any(EvUID == ID_MAIN, axis=1)
            has1 = ak.any(EvUID == det_id, axis=1)

            coinc_mask = has0 & has1

            EvE   = EvE[coinc_mask]
            EvUID = EvUID[coinc_mask]
            #--------------------

            # Create event-level masks for UID == 0 and UID == 1
            mask0 = EvUID == ID_MAIN
            mask1 = EvUID == det_id

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
        ax[det_idx].scatter(sum_EvE_0, sum_EvE_1, color='green', s=5)
        ax[det_idx].set_title(f"{energy_kev} Coincidence between (Ev.UID): {ID_MAIN} and {det_id}, number of coincidence: {n_co}", fontsize=10)
        ax[det_idx].grid(True, alpha=0.5)
        ax[det_idx].set_ylabel(f"E{ID_MAIN}")
        ax[det_idx].set_xlabel(f"E{det_id}")       
        

fig.suptitle(f"Coincidence", fontsize=14, y=0.98)
print("Processing done, witing for plots...")
# --- 5. Final Display ---

plt.show()
