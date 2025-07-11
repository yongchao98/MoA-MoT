import numpy as np

def count_odmr_resonances():
    """
    Calculates the number of ODMR resonances for a hypothetical NV-like center
    in a cubic lattice under an electric field parallel to a cubic edge.

    The logic is based on the Stark effect's dependence on the relative orientation
    of the electric field and the NV center's symmetry axis.
    """

    # 1. In a cubic lattice, assume the three principal axes <100> are the possible
    #    orientations for the NV center's symmetry axis.
    orientations = [np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1])]
    orientation_names = ['[100]', '[010]', '[001]']

    # 2. The electric field is applied along one cubic edge. Let's choose the z-axis.
    e_field_direction = np.array([0, 0, 1])

    # We will use a set to store the labels of unique resonance frequencies.
    # D: zero-field splitting frequency
    # d_perp: splitting magnitude for a transverse E-field
    # d_para_shift: frequency shift for a parallel E-field
    unique_resonances = set()

    print("Analyzing ODMR resonances for an NV-like center in a cubic lattice...")
    print(f"Applied Electric Field is along: {list(e_field_direction)}\n")

    # 3. & 4. Iterate through each possible defect orientation and count resonances.
    for i, orientation in enumerate(orientations):
        # Use the dot product to check for parallel or perpendicular alignment.
        # Since vectors are normalized, dot_product is cos(theta).
        dot_product = np.dot(orientation, e_field_direction)

        print(f"--- Population with orientation {orientation_names[i]} ---")
        if np.isclose(dot_product, 0):
            # Perpendicular case (dot_product is 0)
            print("  - The E-field is PERPENDICULAR to the defect axis.")
            print("  - Effect: The m_s = +/-1 degeneracy is lifted, creating 2 resonances.")
            print("  - Resonances at frequencies: D - d_perp and D + d_perp")
            unique_resonances.add("D - d_perp")
            unique_resonances.add("D + d_perp")
        elif np.isclose(abs(dot_product), 1):
            # Parallel case (dot_product is 1 or -1)
            print("  - The E-field is PARALLEL to the defect axis.")
            print("  - Effect: The m_s = +/-1 states are not split, resulting in 1 resonance.")
            print("  - Resonance at frequency: D + d_para_shift")
            unique_resonances.add("D + d_para_shift")
        print()
    
    # 5. The total number of unique resonances is the size of the set.
    num_resonances = len(unique_resonances)
    
    print("="*40)
    print("Summary:")
    print("The unique resonance frequencies observed from all populations are:")
    # Sort for consistent output
    for resonance in sorted(list(unique_resonances)):
        print(f"  - {resonance}")
    print(f"\nFinal count of distinct ODMR resonances = {num_resonances}")


if __name__ == '__main__':
    count_odmr_resonances()