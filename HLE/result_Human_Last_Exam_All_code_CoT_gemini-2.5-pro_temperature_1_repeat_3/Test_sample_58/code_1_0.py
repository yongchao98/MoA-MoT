import sys

def solve_topological_classification():
    """
    Calculates the topological invariant group for a free fermion system with a defect,
    based on the tenfold way classification.
    """

    # --- Problem Definition ---
    # Symmetries of the 2D free fermion model
    T_squared = -1
    P_squared = -1
    
    # The bulk system is 2-dimensional
    system_spatial_dimension = 2
    
    # The prompt specifies a "point defect (codimension D=1)". This is contradictory
    # for a 2D system. A point defect (0-dimensional) has codimension 2-0=2.
    # A defect with codimension 1 is a line defect (1-dimensional).
    # We will assume "point defect" is the intended type, making the codimension 2.
    defect_dimension = 0

    # --- Step 1: Identify the AZ Symmetry Class ---
    # We map the T and P symmetry properties to one of the 10 AZ classes.
    def get_symmetry_class(tsq, psq):
        if tsq == -1 and psq == -1:
            return "CII"
        elif tsq == 1 and psq == 1:
            return "BDI"
        elif tsq == -1 and psq == 1:
            return "DIII"
        elif tsq == 1 and psq == -1:
            return "CI"
        elif tsq == -1:
            return "AII"
        elif tsq == 1:
            return "AI"
        elif psq == -1:
            return "C"
        elif psq == 1:
            return "D"
        # The remaining are chiral classes A and AIII, which are not determined by T^2/P^2 alone.
        else:
            return "Unknown"

    symmetry_class = get_symmetry_class(T_squared, P_squared)
    if symmetry_class == "Unknown":
        print("Could not determine symmetry class.", file=sys.stderr)
        return

    # --- Step 2: Determine the Defect Codimension ---
    codimension = system_spatial_dimension - defect_dimension

    # --- Step 3: Use the Periodic Table ---
    # The classification of defects of codimension `d` is given by the classification
    # of bulk systems of dimension `d`.
    # Source: Ryu, Schnyder, Furusaki, Ludwig (2010)
    periodic_table = {
        # AZ Class: {dimension: Invariant Group}
        "A":    {0: "Z", 1: "0", 2: "Z", 3: "0", 4: "Z", 5: "0", 6: "Z", 7: "0"},
        "AIII": {0: "0", 1: "Z", 2: "0", 3: "Z", 4: "0", 5: "Z", 6: "0", 7: "Z"},
        "AI":   {0: "Z", 1: "0", 2: "0", 3: "0", 4: "Z", 5: "0", 6: "Z2", 7: "Z2"},
        "BDI":  {0: "Z2", 1: "Z", 2: "0", 3: "0", 4: "0", 5: "Z", 6: "0", 7: "Z2"},
        "D":    {0: "Z2", 1: "Z2", 2: "Z", 3: "0", 4: "0", 5: "0", 6: "Z", 7: "0"},
        "DIII": {0: "0", 1: "Z2", 2: "Z2", 3: "Z", 4: "0", 5: "0", 6: "0", 7: "Z"},
        "AII":  {0: "Z", 1: "0", 2: "Z2", 3: "Z2", 4: "Z", 5: "0", 6: "0", 7: "0"},
        "CII":  {0: "0", 1: "Z", 2: "0", 3: "Z2", 4: "Z2", 5: "Z", 6: "0", 7: "0"},
        "C":    {0: "0", 1: "0", 2: "Z", 3: "0", 4: "Z2", 5: "Z2", 6: "Z", 7: "0"},
        "CI":   {0: "0", 1: "0", 2: "0", 3: "Z", 4: "0", 5: "Z2", 6: "Z2", 7: "Z"},
    }

    # --- Step 4: Find the Invariant Group ---
    try:
        # The table is periodic with period 8 for real classes
        dimension_for_lookup = codimension % 8
        invariant_group = periodic_table[symmetry_class][dimension_for_lookup]
    except KeyError:
        print(f"Could not find entry for Class {symmetry_class} at dimension {dimension_for_lookup}", file=sys.stderr)
        return
        
    # --- Final Output ---
    # The output format explains how the input parameters lead to the final result.
    print("--- Topological Invariant Calculation ---")
    print(f"Input Symmetries: T^2 = {T_squared}, P^2 = {P_squared}")
    print(f"System: 2D with a point defect (dimension={defect_dimension})")
    print("-----------------------------------------")
    print(f"Derived Symmetry Class: {symmetry_class}")
    print(f"Derived Defect Codimension (d): {codimension}")
    print("-----------------------------------------")
    print(f"Topological Invariant Group (from table for class {symmetry_class}, d={codimension}): {invariant_group}")

if __name__ == '__main__':
    solve_topological_classification()