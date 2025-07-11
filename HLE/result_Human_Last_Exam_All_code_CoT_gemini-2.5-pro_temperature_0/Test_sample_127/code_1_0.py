import numpy as np

def solve_odmr_resonances():
    """
    Calculates the number of ODMR resonances for a hypothetical NV center
    in a cubic lattice with an electric field applied along one edge.
    """

    # 1. Define the possible NV center orientations in a cubic lattice (<100> directions)
    # We represent them as unit vectors along X, Y, and Z.
    orientations = {
        "Group 1 (X-axis)": np.array([1, 0, 0]),
        "Group 2 (Y-axis)": np.array([0, 1, 0]),
        "Group 3 (Z-axis)": np.array([0, 0, 1]),
    }

    # 2. Define the electric field direction (parallel to one edge, the X-axis).
    e_field_direction = np.array([1, 0, 0])
    print(f"Electric field is applied along the direction: {e_field_direction}\n")

    # 3. Use a set to store the unique, symbolic resonance frequencies.
    # 'D' represents the zero-field splitting frequency.
    # 'dE' represents the Stark shift.
    distinct_resonances = set()

    print("Analyzing the resonance lines for each NV orientation group:")
    print("="*60)

    # 4. Analyze each orientation group
    for name, nv_axis in orientations.items():
        # The Stark splitting is proportional to the E-field component perpendicular to the NV axis.
        # We can check for parallel vs. perpendicular alignment using the dot product.
        dot_product = np.dot(nv_axis, e_field_direction)

        print(f"NV Group: {name}, Axis: {nv_axis}")

        if dot_product == 1:
            # E-field is PARALLEL to the NV axis.
            # Perpendicular E-field component is zero, so no Stark splitting.
            num_lines = 1
            resonances = ["D"]
            print(f"  - E-field is PARALLEL to the NV axis.")
            print(f"  - This group contributes {num_lines} resonance line at frequency: D")

        elif dot_product == 0:
            # E-field is PERPENDICULAR to the NV axis.
            # This causes a Stark splitting of the m_s = +/-1 levels.
            num_lines = 2
            resonances = ["D - dE", "D + dE"]
            print(f"  - E-field is PERPENDICULAR to the NV axis.")
            print(f"  - This group contributes {num_lines} resonance lines at frequencies: D - dE and D + dE")

        # Add the resulting resonance frequencies to our set of distinct resonances
        for r in resonances:
            distinct_resonances.add(r)

        print("-" * 60)

    # 5. Count the total number of distinct resonance lines
    total_distinct_lines = len(distinct_resonances)

    print("\n--- FINAL RESULT ---")
    print(f"The set of all distinct resonance frequencies is: {sorted(list(distinct_resonances))}")
    print(f"The total number of distinct resonances observed in the ODMR spectrum is: {total_distinct_lines}")

# Run the calculation
solve_odmr_resonances()