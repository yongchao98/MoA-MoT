import sys

def solve_defect_classification():
    """
    Calculates the topological invariant group for a defect in a 2D free fermion system
    with T^2=-1 and P^2=-1 symmetries and a defect of codimension D=1.
    """

    # Step 1: Define system parameters from the problem description
    d_bulk = 2
    T2 = -1
    P2 = -1
    D = 1
    
    print("Step 1: Identify System Parameters")
    print(f"Bulk spatial dimension: d_bulk = {d_bulk}")
    print(f"Time-reversal symmetry: T^2 = {T2}")
    print(f"Particle-hole symmetry: P^2 = {P2}")
    print(f"Defect codimension: D = {D}")
    print("-" * 30)

    # Step 2: Determine the symmetry class
    # The combination T^2=-1, P^2=-1 defines the AZ class CII.
    # We map the 10 real AZ classes to an index s = 0..7.
    class_to_s = {
        "AI": 0, "BDI": 1, "D": 2, "DIII": 3, 
        "AII": 4, "CII": 5, "C": 6, "CI": 7
    }
    s_to_class = {v: k for k, v in class_to_s.items()}
    
    original_class = "CII"
    s_original = class_to_s[original_class]
    
    print("Step 2: Determine the Symmetry Class")
    print(f"The symmetries T^2={T2} and P^2={P2} correspond to the Altland-Zirnbauer class: {original_class}")
    print(f"This class has the index s = {s_original} in the periodic table.")
    print("-" * 30)

    # Step 3: Apply the defect classification formula
    # The classification group for a defect of dimension d_defect in class s_original
    # is equivalent to the classification of a 0D system in class s_new.
    # d_defect = d_bulk - D
    # s_new = (s_original - d_defect) mod 8
    
    d_defect = d_bulk - D
    s_new = (s_original - d_defect) % 8
    new_class = s_to_class[s_new]
    
    print("Step 3: Apply Defect Classification Theory")
    print("The problem is mapped to finding the invariant for a 0D system in a new symmetry class.")
    print(f"First, we calculate the dimension of the defect:")
    print(f"  d_defect = d_bulk - D = {d_bulk} - {D} = {d_defect}")
    print("Next, we calculate the index 's' of the new symmetry class:")
    print(f"  s_new = (s_original - d_defect) mod 8 = ({s_original} - {d_defect}) mod 8 = {s_new}")
    print(f"The class with index s_new = {s_new} is: {new_class}")
    print("-" * 30)

    # Step 4: Look up the result in the periodic table
    # We need the classification for a 0D system (d=0) in the new class.
    # The table stores invariants for d=0 to d=7.
    periodic_table = {
        # s: (d=0, d=1, d=2, d=3, d=4, d=5, d=6, d=7)
        "AI":   ("Z", "Z2", "Z2", "0", "Z", "0", "0", "0"),
        "BDI":  ("Z", "Z", "Z2", "Z2", "0", "Z", "0", "0"),
        "D":    ("Z2", "Z", "Z", "Z2", "Z2", "0", "Z", "0"),
        "DIII": ("Z2", "Z2", "Z", "Z", "Z2", "Z2", "0", "Z"),
        "AII":  ("Z", "Z2", "Z2", "Z", "Z", "Z2", "Z2", "0"),
        "CII":  ("0", "Z", "Z2", "Z2", "Z", "Z", "Z2", "Z2"),
        "C":    ("0", "0", "Z", "Z2", "Z2", "Z", "Z", "Z2"),
        "CI":   ("0", "0", "0", "Z", "Z2", "Z2", "Z", "Z")
    }

    # The classification is the d=0 entry for the new_class
    final_invariant_group = periodic_table[new_class][0]

    print("Step 4: Find the Invariant in the Periodic Table")
    print(f"We need the classification for a d=0 system in class {new_class}.")
    print(f"Looking up periodic_table['{new_class}'][0] gives the result.")
    print("\nFinal Answer:")
    print(f"The topological invariant group is: {final_invariant_group}")

if __name__ == "__main__":
    solve_defect_classification()