import sys

def solve_tenfold_way_defect():
    """
    Calculates the topological invariant group for a defect in a free-fermion system
    based on the tenfold way classification.
    """
    # Step 1: Define the parameters from the problem description.
    # d: spatial dimension of the bulk system
    # D: codimension of the defect
    # T_sq: Value of T^2 (Time-Reversal symmetry)
    # P_sq: Value of P^2 (Particle-Hole symmetry)
    d = 2
    D = 1
    T_sq = -1
    P_sq = -1

    print("--- Problem Parameters ---")
    print(f"Bulk spatial dimension: d = {d}")
    print(f"Defect codimension: D = {D}")
    print(f"Symmetries: T^2 = {T_sq}, P^2 = {P_sq}\n")

    # Step 2: Determine the Altland-Zirnbauer (AZ) symmetry class.
    # The 8 "real" classes are indexed by q from 0 to 7.
    # We create a mapping from (T^2, P^2) to the class name and index q.
    real_az_classes = {
        (1, 0):   ("AI", 0),
        (1, 1):   ("BDI", 1),
        (0, 1):   ("D", 2),
        (-1, 1):  ("DIII", 3),
        (-1, 0):  ("AII", 4),
        (-1, -1): ("CII", 5),
        (0, -1):  ("C", 6),
        (1, -1):  ("CI", 7),
    }

    key = (T_sq, P_sq)
    if key not in real_az_classes:
        print(f"Error: Symmetry combination (T^2={T_sq}, P^2={P_sq}) is not a standard real AZ class.")
        sys.exit(1)

    class_name, q = real_az_classes[key]
    print("--- Step 1: Identify Symmetry Class ---")
    print(f"The system belongs to the Altland-Zirnbauer class '{class_name}'.")
    print(f"This class corresponds to the real class index q = {q}.\n")

    # Step 3: Calculate the effective dimension 's'.
    s = d - D
    print("--- Step 2: Calculate Effective Dimension ---")
    print("The topological invariant of a defect is determined by its effective dimension 's'.")
    print(f"The calculation is: s = d - D")
    print(f"Numerically: s = {d} - {D} = {s}\n")

    # Step 4: Apply Bott periodicity to find the equivalent 0D class.
    # The classification in dimension 's' for class 'q' is given by the
    # classification in dimension 0 for a class with index k = (q - s) mod 8.
    k = (q - s) % 8
    print("--- Step 3: Apply Bott Periodicity ---")
    print("The classification is periodic. We can map the problem to an equivalent 0D problem.")
    print("The new class index 'k' is calculated as: k = (q - s) mod 8")
    print(f"Numerically: k = ({q} - {s}) mod 8 = {k}\n")

    # Step 5: Look up the 0D classification for class k.
    # This is the pi_0 group for the classifying space R_k.
    # k: (Class Name, Invariant Group)
    d0_classification = {
        0: ("AI", "Z"),
        1: ("BDI", "Z_2"),
        2: ("D", "Z_2"),
        3: ("DIII", "0"),
        4: ("AII", "2Z"),
        5: ("CII", "0"),
        6: ("C", "0"),
        7: ("CI", "0"),
    }
    
    k_class_name, final_group = d0_classification[k]
    print("--- Step 4: Find the 0D Invariant Group ---")
    print(f"The classification is given by the invariant of class k={k} (which is {k_class_name}) in d=0.")
    print(f"The invariant group for class {k_class_name} in 0 dimensions is: {final_group}\n")

    # Conclusion
    print("--- Final Result ---")
    print(f"The group of the topological invariant for this system is {final_group}.")
    print("(Z: integers, Z_2: cyclic group of order 2, 2Z: even integers, 0: trivial group)")

solve_tenfold_way_defect()