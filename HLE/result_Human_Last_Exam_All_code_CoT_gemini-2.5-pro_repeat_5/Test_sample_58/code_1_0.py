def solve_topological_classification():
    """
    Determines the topological invariant group for a 2D free fermion system
    with specified symmetries and a codimension 1 defect.
    """

    # --- Problem Parameters ---
    # The system is a 2D free fermion model.
    spatial_dimension = 2
    # It has time-reversal symmetry T with T^2 = -1.
    time_reversal_sq = -1
    # It has particle-hole symmetry P with P^2 = -1.
    particle_hole_sq = -1
    # The system contains a defect of codimension D = 1.
    defect_codimension = 1

    # --- Step 1: Identify the Symmetry Class ---
    # The 10-fold way classifies systems based on T, P, and Chiral symmetries.
    # The combination T^2 = -1 and P^2 = -1 corresponds to the AZ class CII.
    symmetry_class = "CII"

    # --- Step 2: Apply Defect Classification Rule ---
    # The classification of a defect of codimension D is equivalent to the
    # classification of a bulk system with dimension d' = D.
    effective_dimension = defect_codimension

    # --- Step 3: Look up in the Periodic Table ---
    # We use a data structure to represent the periodic table for the 8 real AZ classes.
    # The table gives the classification (Z, Z2, or 0) for dimensions d = 0 to 7.
    # 0: Trivial group {0}
    # Z: Integers
    # Z2: Integers modulo 2
    periodic_table = {
        # class: [d=0, d=1, d=2, d=3, d=4, d=5, d=6, d=7]
        "BDI":   ["Z", "Z2", "Z2", "0", "Z", "0", "0", "0"],
        "D":     ["Z2", "Z2", "0", "Z", "0", "0", "0", "Z"],
        "DIII":  ["Z2", "0", "Z", "0", "0", "0", "Z", "Z2"],
        "AII":   ["0", "Z", "0", "0", "0", "Z", "Z2", "Z2"],
        "CII":   ["0", "0", "0", "Z", "Z2", "Z2", "0", "Z"],
        "C":     ["0", "0", "Z", "Z2", "Z2", "0", "Z", "0"],
        "CI":    ["0", "Z", "Z2", "Z2", "0", "Z", "0", "0"],
        "AI":    ["Z", "0", "0", "0", "Z", "Z2", "Z2", "0"],
    }
    
    # Retrieve the result from the table
    classification_group = periodic_table[symmetry_class][effective_dimension]

    # --- Step 4: Print the Final Result ---
    print("--- Topological Defect Classification ---")
    print(f"1. System Parameters:")
    print(f"   - Spatial dimension d = {spatial_dimension}")
    print(f"   - Time-reversal property T^2 = {time_reversal_sq}")
    print(f"   - Particle-hole property P^2 = {particle_hole_sq}")
    print(f"   - Defect codimension D = {defect_codimension}")
    print("\n2. Analysis:")
    print(f"   - The specified symmetries place the system in AZ class '{symmetry_class}'.")
    print(f"   - The classification for a defect of codimension D={defect_codimension} is equivalent to a bulk system of dimension d'={effective_dimension}.")
    print(f"   - Looking up class '{symmetry_class}' at dimension d'={effective_dimension} in the periodic table gives the group '{classification_group}'.")
    print("\n3. Conclusion:")
    print(f"The group of the topological invariant is: {classification_group}")


solve_topological_classification()