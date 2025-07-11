import sys

def solve_tenfold_way_defect():
    """
    Calculates the topological invariant group for a specified fermionic system with a defect.
    """
    # Step 1: Define system properties from the user's request.
    T_sq = -1
    P_sq = -1
    codimension_D = 1

    # In the tenfold way, the combination of T and P symmetries determines the class.
    # T^2=-1, P^2=-1 implies Chiral symmetry S^2=1. This is class CII.
    # We assign an index 'q' to each of the 8 real AZ classes.
    # (T^2, P^2) -> Class Name (q)
    # (+1, +1) -> BDI (1)
    # (-1, +1) -> DIII (3)
    # (-1, -1) -> CII (5) <--- This is our case.
    # (+1, -1) -> CI (7)
    # Others have T or P being 0.
    class_name = "CII"
    q_index = 5

    # Step 2: Apply the formula for defect classification.
    # The group K is the (D-1)-th homotopy group of the classifying space R_q.
    # K = pi_{D-1}(R_q)
    homotopy_k = codimension_D - 1

    # Step 3: Relate to the periodic table using Bott periodicity.
    # The group pi_k(R_q) is isomorphic to the bulk classification of a 0D system
    # in the class C_{q-k}. We find the table entry for the effective index (q-k) mod 8.
    effective_index = (q_index - homotopy_k) % 8

    # The periodic table for bulk 0D systems (d=0), indexed by `q-d`.
    # Index (q-0) | Class | Group
    # ---------------------------
    #      0       |  AI   |   Z
    #      1       |  BDI  |  Z_2
    #      2       |   D   |  Z_2
    #      3       | DIII  |   0
    #      4       |  AII  |   Z
    #      5       |  CII  |   0   <--- This is our case.
    #      6       |   C   |   0
    #      7       |   CI  |   0
    classification_table_d0 = {
        0: "Z", 1: "Z_2", 2: "Z_2", 3: "0",
        4: "Z", 5: "0", 6: "0", 7: "0"
    }
    result_group = classification_table_d0[effective_index]

    # Step 4: Print the full derivation and final result.
    print("Derivation of the Topological Invariant Group:")
    print("-" * 45)
    print(f"1. The system has symmetries T^2 = {T_sq} and P^2 = {P_sq}.")
    print(f"   This corresponds to the Altland-Zirnbauer (AZ) class '{class_name}', which has the index q = {q_index}.")
    print(f"2. The defect is specified to have codimension D = {codimension_D}.")
    print(f"3. The classification group 'K' is given by the homotopy group K = \u03c0_{{D-1}}(R_q).")
    print(f"   Substituting the given values: K = \u03c0_{{{codimension_D}-1}}(R_{q_index}) = \u03c0_{homotopy_k}(R_{q_index}).")
    print(f"4. Using the periodicity theorem, this is equivalent to the bulk classification of a 0D system")
    print(f"   for the class with an effective index of (q - k) mod 8 = ({q_index} - {homotopy_k}) mod 8 = {effective_index}.")
    print(f"5. The periodic table shows the group for this effective index is '{result_group}'.")
    print("-" * 45)
    
    # Final equation as requested by the prompt.
    print("\nFinal Equation and Answer:")
    print(f"Invariant Group = \u03c0_({codimension_D} - 1)(R_{q_index}) = \u03c0_{homotopy_k}(R_{q_index}) = {result_group}")

solve_tenfold_way_defect()