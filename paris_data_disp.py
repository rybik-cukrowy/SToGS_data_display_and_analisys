import uproot
import numpy as np
import matplotlib.pyplot as plt
import awkward as ak
import os

# CONSTANTS

N_DET = 52

SOURCE_ENERGIES_KEV = [200, 2000, 10000]

TREE_KEYS = ["SToGS;1", "SToGS;2"]

BIN_MAP = { # used for histogram displays
    200:  np.linspace(0, 0.2, 20),
    2000: np.linspace(0, 2, 300),
    10000: np.linspace(0, 11, 300),
}

det_groups = [
        range(0, 50,2), # parzyste id
        range(1, 51, 2), # nieparzyste id
    ]

fig, ax = plt.subplots(len(SOURCE_ENERGIES_KEV), len(det_groups))
fig.subplots_adjust(hspace=0.4)


# DATA

for row_idx, energy_kev in enumerate(SOURCE_ENERGIES_KEV): # one plot (two detectors) per iteration
    # FILE NAMING DEFINITION
    filename = f"paris1//paris1_{energy_kev}.root"
    print(f"Processing file: {filename}...")
    

    if not os.path.exists(filename):
        print(f"WARNING: File not found: {filename}. Skipping this energy.")
        continue
        
    for tree_key in TREE_KEYS:

        # part that is acrually doing something
        with uproot.open(filename) as file:
            for name, classname in file.classnames().items():
                if classname == "TTree":
                    print(name)
            tree = file[tree_key]
             

            EvE = tree["Ev.E"].array()  # energies of hits
            EvUID = tree["Ev.UID"].array()  # number of hits per event

        

        per_detector_energies = []

        for det in range(N_DET):
            # 1. Select hits from this detector (ragged-safe)
            det_hits = EvE[EvUID == det]

            # 2. Sum hits per event (this is the crucial step)
            det_event_energy = ak.sum(det_hits, axis=1)

            # 3. Keep only events where this detector actually saw energy
            det_event_energy = det_event_energy[det_event_energy > 0]

            # 4. Convert to NumPy: this is your "box of energies"
            per_detector_energies.append(
                ak.to_numpy(det_event_energy)
            )
        
    # Plotting the histograms: 13 separate lines per batch
    
    for group_idx, group in enumerate(det_groups):
        for detector in group:
            ax[row_idx, group_idx].hist(
                per_detector_energies[detector],
                bins=BIN_MAP[energy_kev],
                histtype="step",
                label=f"Det {detector}"
            )
        ax[row_idx, group_idx].set_title(
            f"{energy_kev/1000:.1f} MeV source, group: {group_idx}"
        )
        ax[row_idx, group_idx].set_yscale("log")
        ax[row_idx, group_idx].grid("True")
        #ax[row_idx, group_idx].legend()
        



fig.suptitle(f"PARIS array", fontsize=14, y=0.98)
print("Processing done, witing for plots...")
# --- 5. Final Display ---

plt.show()