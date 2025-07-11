import numpy as np

def solve_odmr_resonances():
    """
    Calculates the number of ODMR resonances for a hypothetical NV center in a cubic lattice
    under an electric field parallel to a lattice edge.
    """
    # 1. Define the possible NV axis orientations for a cubic lattice.
    # These are the three principal axes <100>.
    nv_orientations = {
        "Group 1": np.array([1, 0, 0]),
        "Group 2": np.array([0, 1, 0]),
        "Group 3": np.array([0, 0, 1])
    }

    # 2. Define the electric field direction, parallel to a lattice edge.
    E_field_direction = np.array([1, 0, 0])

    print("Analyzing ODMR resonances for a hypothetical NV center in a cubic lattice.")
    print(f"Electric field is applied along the {E_field_direction} direction.\n")

    # 3. Use a set to store the symbolic names of the unique resonance frequencies.
    # Let 'D' be the zero-field splitting frequency.
    # Let 'd' be the Stark splitting magnitude.
    unique_resonances = set()
    num_peaks_per_group = {}

    for group_name, nv_axis in nv_orientations.items():
        # Calculate the dot product to check for alignment.
        # cos(theta) = (E · n) / (|E| |n|)
        # Since vectors are normalized, dot product is cos(theta).
        dot_product = np.dot(E_field_direction, nv_axis)

        print(f"--- {group_name} (NV axis || {nv_axis}) ---")
        if np.isclose(dot_product, 1.0):
            # E-field is parallel to the NV axis.
            # Perpendicular component is zero, so no splitting.
            num_peaks = 1
            unique_resonances.add("D")
            num_peaks_per_group[group_name] = num_peaks
            print(f"E-field is PARALLEL to the NV axis. This group contributes {num_peaks} resonance.")
        elif np.isclose(dot_product, 0.0):
            # E-field is perpendicular to the NV axis.
            # Perpendicular component is maximum, causing splitting.
            num_peaks = 2
            unique_resonances.add("D - d")
            unique_resonances.add("D + d")
            num_peaks_per_group[group_name] = num_peaks
            print(f"E-field is PERPENDICULAR to the NV axis. This group contributes {num_peaks} resonances.")
        else:
            # This case doesn't happen for the chosen E-field and orientations, but is included for completeness.
            num_peaks = 2
            # Splitting would depend on the angle.
            unique_resonances.add("D - d'")
            unique_resonances.add("D + d'")
            num_peaks_per_group[group_name] = num_peaks
            print(f"E-field is at an angle to the NV axis. This group contributes {num_peaks} resonances.")
        print("-" * 45)

    # 4. Summarize and present the final count.
    print("\n--- SUMMARY ---")
    print(f"The set of unique resonance frequencies is: {sorted(list(unique_resonances))}")

    # The final count is the number of unique resonances found.
    # It can be thought of as the contribution from the parallel group plus the contribution from the perpendicular groups.
    num_parallel_peaks = 1 # From Group 1
    num_perpendicular_peaks = 2 # From Group 2 and 3, which are degenerate
    total_distinct_peaks = len(unique_resonances)

    print("\nThe final equation for the total number of distinct resonance lines is:")
    print(f"Total = (peaks from parallel group) + (peaks from perpendicular groups)")
    print(f"Total = {num_parallel_peaks} + {num_perpendicular_peaks} = {total_distinct_peaks}")

    print(f"\nTherefore, you will see a total of {total_distinct_peaks} resonances in the ODMR spectrum.")

    return total_distinct_peaks

# Execute the function
solve_odmr_resonances()