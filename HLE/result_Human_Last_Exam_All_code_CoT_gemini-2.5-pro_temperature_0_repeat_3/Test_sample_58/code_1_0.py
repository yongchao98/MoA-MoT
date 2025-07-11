import sys

def solve_fermion_defect_classification():
    """
    Calculates the topological invariant group for a defect in a 2D free fermion system.
    """
    # Step 1: Define system parameters from the user's query.
    d = 2   # Spatial dimension
    T2 = -1 # Time-reversal symmetry squared
    P2 = -1 # Particle-hole symmetry squared

    # The prompt states "point defect (codimension D=1)". This is ambiguous.
    # A point defect in 2D has codimension q=2. A line defect has q=1.
    # We will proceed with the explicitly stated codimension q=1.
    q = 1

    # Step 2: Determine the symmetry class and its index 's'.
    # We use a dictionary to map (T^2, P^2) to the class name and index.
    # A value of 0 for T2 or P2 indicates the absence of that symmetry.
    real_classes = {
        # (T^2, P^2): (class_name, s_index)
        (1, 0): ('AI', 0), (1, 1): ('BDI', 1), (0, 1): ('D', 2),
        (-1, 1): ('DIII', 3), (-1, 0): ('AII', 4), (-1, -1): ('CII', 5),
        (0, -1): ('C', 6), (1, -1): ('CI', 7),
    }

    # Find the class based on T2 and P2.
    # For real classes, if T and P are present, Chiral symmetry S=TP is also present.
    key = (T2, P2)
    if key not in real_classes:
        print(f"Error: Symmetry combination (T^2={T2}, P^2={P2}) is not a standard real class.", file=sys.stderr)
        return
    
    class_name, s = real_classes[key]

    # Step 3 & 4: Apply K-theory formula and Bott periodicity.
    # The classification group is π_{q-1}(R_{d-s}), which simplifies to π_0(R_{d-s+q-1}).
    # We calculate the final index n = (d - s + q - 1) mod 8.
    index = (d - s + q - 1) % 8

    # The sequence of classification groups π_0(R_n) for n = 0, 1, ..., 7.
    pi0_R_n = {
        0: 'Z', 1: 'Z_2', 2: 'Z_2', 3: '0',
        4: 'Z', 5: '0', 6: '0', 7: '0',
    }
    
    group = pi0_R_n[index]

    # Step 5: Print the step-by-step derivation and the final result.
    print("--- Calculation of the Topological Invariant Group ---")
    print(f"1. System Identification:")
    print(f"   - Spatial Dimension (d): {d}")
    print(f"   - Symmetries: T^2 = {T2}, P^2 = {P2}")
    print(f"   - Defect Codimension (q): {q} (Interpreted as a line defect)")
    print(f"   - This corresponds to symmetry class '{class_name}' with index s = {s}.")
    
    print("\n2. Classification Formula:")
    print(f"   - The group is given by π_{{q-1}}(R_{{d-s}}), which via Bott periodicity is π_0(R_{{d-s+q-1}}).")
    print(f"   - We calculate the final index 'n' for the table lookup.")
    
    print("\n3. Final Equation:")
    # The user requested to output each number in the final equation.
    print(f"   n = (d - s + q - 1) mod 8")
    print(f"   n = ({d} - {s} + {q} - 1) mod 8 = {index}")

    print("\n4. Result:")
    print(f"   The group is π_0(R_{{{index}}}), which from the periodic table is '{group}'.")
    print("------------------------------------------------------")

    # For completeness, check the other interpretation of the prompt.
    q_alt = 2 # Codimension for a point defect in 2D
    index_alt = (d - s + q_alt - 1) % 8
    group_alt = pi0_R_n[index_alt]
    print("\nNote: If 'point defect' (q=2) was intended, the result is the same:")
    print(f"   n = ({d} - {s} + {q_alt} - 1) mod 8 = {index_alt}, giving the group '{group_alt}'.")

    # The final answer is the trivial group, denoted by 0.
    # print(f"\n<<<{group}>>>")

if __name__ == '__main__':
    solve_fermion_defect_classification()