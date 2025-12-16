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
            # arraye do energi na jeden symulacje
            T1_avg = []
            T0_avg = []

            print(f"Processing tree: {tree_key}...")

            tree = uproot.open(filename)[tree_key]
            EvE   = tree["Ev.E"].array()
            EvUID = tree["Ev.UID"].array()
            EvT   = tree["Ev.T"].array()

            # coincidence mask
            has0 = ak.any(EvUID == 0, axis=1)
            has1 = ak.any(EvUID == 1, axis=1)
            coinc_mask = has0 & has1

            EvE   = EvE[coinc_mask]
            EvT   = EvT[coinc_mask]
            EvUID = EvUID[coinc_mask]

            n_co = len(EvE)

            # detector masks
            mask0 = EvUID == 0
            mask1 = EvUID == 1

            # split by detector
            EvE_0 = EvE[mask0]
            EvE_1 = EvE[mask1]
            EvT_0 = EvT[mask0]
            EvT_1 = EvT[mask1]

            # energy sums (if needed)
            sum_EvE_0 = ak.sum(EvE_0, axis=1)
            sum_EvE_1 = ak.sum(EvE_1, axis=1)

            # energy-weighted average times
            T0_avg = ak.sum(EvT_0 * EvE_0, axis=1) / ak.sum(EvE_0, axis=1)
            T1_avg = ak.sum(EvT_1 * EvE_1, axis=1) / ak.sum(EvE_1, axis=1)
          
        

        # Plotting the histograms
        T0_np = ak.to_numpy(T0_avg)
        T1_np = ak.to_numpy(T1_avg)

        delta_T = T1_np-T0_np

        bins = np.linspace(-2, 2, 500) # wartości dopasowanie bo tak

        ax[row_idx, orient_idx].hist(
            delta_T,
            bins=bins,
            histtype="step",
            label=f"{energy_kev} keV",
            color="lime"
        )

        ax[row_idx, orient_idx].set_title(
            f"{energy_kev/1000} MeV, orientation {orientation}, coincidences: {n_co}",
            fontsize=10
        )
        ax[row_idx,orient_idx].grid(True, alpha=0.5)
        # ax[row_idx,orient_idx].set_yscale("log")
        ax[row_idx,orient_idx].set_ylabel("Counts")
        ax[row_idx,orient_idx].set_xlabel("delta T [ns]")
       
        


fig.suptitle(f"Detection time difference between different detector relative placement", fontsize=14, y=0.98)
print("Processing done, witing for plots...")
# --- 5. Final Display ---
plt.show()