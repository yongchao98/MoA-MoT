import numpy as np

def analyze_nv_in_cubic_lattice():
    """
    Calculates the number of ODMR resonances for a hypothetical NV center
    in a cubic lattice with an E-field along a lattice edge.
    """
    # 1. Define the system
    # In a cubic lattice, the three possible NV center orientations are along the
    # principal axes [100], [010], and [001]. We represent them as unit vectors.
    nv_orientations = {
        "along X-axis": np.array([1, 0, 0]),
        "along Y-axis": np.array([0, 1, 0]),
        "along Z-axis": np.array([0, 0, 1]),
    }

    # The electric field is applied parallel to one of the cubic lattice edges.
    # Without loss of generality, we choose the Z-axis.
    e_field_direction = np.array([0, 0, 1])

    print("Analyzing a hypothetical NV center in a cubic lattice.")
    print(f"Electric field direction is set to: {e_field_direction}\n")

    # This set will store the symbolic representations of the unique resonance lines.
    # Using a set automatically handles duplicates from equivalent orientations.
    unique_resonance_formulas = set()
    num_parallel_resonances = 0
    num_perpendicular_resonances = 0


    # 2. Analyze the effect of the E-field on each orientation
    for name, orientation_vec in nv_orientations.items():
        # Check if the E-field is parallel or perpendicular to the NV axis.
        # A dot product of 1 (or -1) means parallel; 0 means perpendicular.
        dot_product = np.abs(np.dot(orientation_vec, e_field_direction))

        if np.isclose(dot_product, 1.0):
            # Case 1: E-field is PARALLEL to the NV axis.
            # This does NOT split the m_s = +/-1 degeneracy.
            # Result: 1 resonance line.
            print(f"For NV {name}: E-field is PARALLEL. Result: 1 resonance.")
            # Symbolic formula for the resonance frequency (D0 is zero-field splitting)
            formula = "f = D0 + delta_parallel"
            unique_resonance_formulas.add(formula)
            num_parallel_resonances = 1 # There is only one resonance in this case

        elif np.isclose(dot_product, 0.0):
            # Case 2: E-field is PERPENDICULAR to the NV axis.
            # This splits the m_s = +/-1 degeneracy.
            # Result: 2 resonance lines.
            print(f"For NV {name}: E-field is PERPENDICULAR. Result: 2 resonances.")
            # Symbolic formulas for the two split resonance frequencies
            formula1 = "f = D0 + delta_perp"
            formula2 = "f = D0 - delta_perp"
            unique_resonance_formulas.add(formula1)
            unique_resonance_formulas.add(formula2)
            # The two perpendicular orientations are equivalent and will produce
            # the same two lines. We only need to count these two lines once.
            num_perpendicular_resonances = 2

    # 3. Count the total number of distinct resonances
    total_resonances = len(unique_resonance_formulas)

    print("\n--- Summary ---")
    print("The unique resonance frequencies observed are:")
    for formula in sorted(list(unique_resonance_formulas)):
        print(f"  - {formula}")

    print("\nFinal Calculation:")
    print("We have one group of NV centers with the E-field parallel to their axis, and two groups with the field perpendicular.")
    print("The two perpendicular groups are equivalent and produce the same pair of resonance lines.")
    print("Therefore, the total number of distinct resonances is the sum from each unique case.")
    print(f"Equation: {num_parallel_resonances} (from parallel case) + {num_perpendicular_resonances} (from perpendicular case) = {total_resonances}")
    print(f"\nTotal distinct ODMR resonances = {total_resonances}")

if __name__ == '__main__':
    analyze_nv_in_cubic_lattice()