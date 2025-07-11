import sys

def solve_topological_classification():
    """
    Solves for the topological invariant group based on the tenfold way classification
    for a system with a defect.
    """
    # --- Problem Parameters ---
    # Spatial dimension of the system
    d = 2
    # Properties of the symmetries (T^2 and P^2)
    T_squared = -1
    P_squared = -1
    # Codimension of the topological defect
    D = 1

    # --- Data Structures for Tenfold Way ---
    # A mapping from (T^2, P^2) values to the Altland-Zirnbauer (AZ) class name.
    # A value of 0 would represent an absent symmetry.
    symmetry_map = {
        (1, 1): "BDI",
        (-1, 1): "DIII",
        (1, -1): "CI",
        (-1, -1): "CII",
        (1, 0): "AI",
        (-1, 0): "AII",
        (0, 1): "D",
        (0, -1): "C",
    }

    # The periodic table for real AZ classes.
    # The value is a list of classification groups for dimensions d = 0 to 7.
    periodic_table_real = {
        # d=0,  d=1,  d=2,   d=3,   d=4,   d=5,  d=6,   d=7
        "AI":   ["Z",  "0",   "0",   "0",   "Z",   "0",  "Z_2", "Z_2"],
        "BDI":  ["Z_2","Z",   "0",   "0",   "0",   "Z",  "0",   "Z_2"],
        "D":    ["Z_2","Z_2", "Z",   "0",   "0",   "0",  "Z",   "0"  ],
        "DIII": ["0",  "Z_2", "Z_2", "Z",   "0",   "0",  "0",   "Z"  ],
        "AII":  ["Z",  "0",   "Z_2", "Z_2", "Z",   "0",  "0",   "0"  ],
        "CII":  ["0",  "Z",   "0",   "Z_2", "Z_2", "Z",  "0",   "0"  ],
        "C":    ["0",  "0",   "Z",   "0",   "Z_2", "Z_2","Z",   "0"  ],
        "CI":   ["0",  "0",   "0",   "Z",   "0",   "Z_2","Z_2", "Z"  ],
    }

    # --- Step-by-Step Solution ---
    print("Solving for the topological invariant group:")
    print(f"System Details:")
    print(f"  - Spatial dimension d = {d}")
    print(f"  - Symmetries: T^2 = {T_squared}, P^2 = {P_squared}")
    print(f"  - Defect codimension D = {D}\n")

    # Step 1: Determine the symmetry class
    class_key = (T_squared, P_squared)
    az_class = symmetry_map.get(class_key)
    if not az_class:
        print("Error: Could not determine the symmetry class for the given inputs.")
        sys.exit(1)
    
    print(f"Step 1: Determine the symmetry class.")
    print(f"The symmetries T^2 = {T_squared} and P^2 = {P_squared} correspond to the Altland-Zirnbauer class '{az_class}'.\n")

    # Step 2: Calculate the effective dimension for the defect classification
    d_eff = d - D
    print(f"Step 2: Apply the bulk-defect correspondence.")
    print(f"The classification for a defect is equivalent to the bulk classification in an effective dimension d' = d - D.")
    print(f"The effective dimension is calculated as: d' = {d} - {D} = {d_eff}\n")

    # Step 3: Look up the result in the periodic table
    # Use Bott periodicity (period of 8 for real classes) if d_eff is negative
    d_lookup = d_eff % 8
    
    classification_group = periodic_table_real[az_class][d_lookup]
    
    print(f"Step 3: Find the classification in the periodic table.")
    print(f"Looking up the entry for class '{az_class}' in dimension d' = {d_eff}...")
    print("--------------------------------------------------")
    print(f"The resulting group of the topological invariant is: {classification_group}")
    print("--------------------------------------------------")

if __name__ == '__main__':
    solve_topological_classification()