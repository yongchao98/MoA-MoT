import numpy as np

def calculate_odmr_resonances():
    """
    Calculates the number of ODMR resonances for a hypothetical NV center
    in a cubic lattice under an electric field along a crystal axis.
    """
    # In a cubic lattice, the NV center's symmetry axis can align with the
    # three principal axes <100>.
    # Let's represent them as unit vectors.
    nv_orientations = {
        "x-axis": np.array([1, 0, 0]),
        "y-axis": np.array([0, 1, 0]),
        "z-axis": np.array([0, 0, 1]),
    }

    # The electric field is applied parallel to one of the cubic lattice edges.
    # Let's choose the x-axis without loss of generality.
    E_field_direction = np.array([1, 0, 0])
    print("Electric field is applied along the x-axis.\n")

    # We use a set to store the unique resonance splittings.
    # 0 represents the unsplit resonance.
    # A non-zero value represents a symmetric split into two resonances.
    unique_splittings = set()

    # The final equation components
    unsplit_count = 0
    split_group_count = 0

    print("Analyzing each possible NV orientation:")
    for name, nv_axis in nv_orientations.items():
        # The splitting of ODMR resonance is proportional to the magnitude of the
        # electric field component perpendicular to the NV axis.
        # E_perp = E - (E . nv_axis) * nv_axis
        E_parallel_component = np.dot(E_field_direction, nv_axis)
        E_perp_magnitude = np.linalg.norm(E_field_direction - E_parallel_component * nv_axis)

        # Use a small tolerance for floating point comparison
        if np.isclose(E_perp_magnitude, 0):
            num_resonances = 1
            splitting_type = "unsplit"
            unique_splittings.add(0)
            unsplit_count += 1
        else:
            num_resonances = 2
            splitting_type = "split"
            # The magnitude of splitting is proportional to E_perp_magnitude
            # We add it to the set to see if this splitting is unique
            unique_splittings.add(round(E_perp_magnitude, 5))
            # We count groups of split peaks, not individual peaks here.
            # The set handles the uniqueness.

        print(f"- For NVs oriented along the {name}:")
        print(f"  - The E-field is {'parallel' if num_resonances == 1 else 'perpendicular'} to the NV axis.")
        print(f"  - This results in {num_resonances} resonance line(s) ({splitting_type}).")

    # The number of split groups is the number of non-zero elements in the set
    split_group_count = len(unique_splittings) - (1 if 0 in unique_splittings else 0)
    
    # Each split group contributes 2 lines
    split_lines_count = split_group_count * 2

    # The number of unsplit lines is 1 if there is a zero in the set
    unsplit_lines_count = 1 if 0 in unique_splittings else 0

    total_resonances = unsplit_lines_count + split_lines_count

    print("\n---")
    print("Summary of ODMR lines:")
    print(f"Number of unsplit resonance lines: {unsplit_lines_count}")
    print(f"Number of split resonance lines: {split_lines_count} (from {split_group_count} unique group(s) of split NVs)")
    print("\nFinal Equation:")
    print(f"{unsplit_lines_count} (unsplit) + {split_lines_count} (split) = {total_resonances} total resonances")
    print(f"\nTherefore, you will see a total of {total_resonances} resonances in the ODMR spectrum.")
    return total_resonances

# Execute the function
calculate_odmr_resonances()