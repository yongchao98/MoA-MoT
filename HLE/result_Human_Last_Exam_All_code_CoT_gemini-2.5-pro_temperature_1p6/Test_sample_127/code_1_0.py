import numpy as np

def calculate_odmr_resonances():
    """
    Calculates the number of ODMR resonances for a hypothetical NV center
    in a cubic lattice under an electric field parallel to a cube edge.
    """

    # 1. Define the possible orientations of the NV center in a cubic lattice.
    # These are along the principal axes <100>.
    nv_orientations = [
        np.array([1, 0, 0]),  # NV axis along x
        np.array([0, 1, 0]),  # NV axis along y
        np.array([0, 0, 1]),  # NV axis along z
    ]

    # 2. Define the applied electric field direction, parallel to a lattice edge.
    e_field_direction = np.array([0, 0, 1]) # Along z-axis

    # 3. Categorize orientations and count resonances.
    # Resonances are represented as strings to handle uniqueness.
    # For a transverse field, splitting gives two peaks: 'D+d' and 'D-d'
    # For a longitudinal field, no splitting gives one peak: 'D_prime'
    resonance_frequencies = set()
    parallel_orientations_count = 0
    perpendicular_orientations_count = 0
    
    for nv_axis in nv_orientations:
        # Use dot product to check for parallel or perpendicular alignment.
        # np.isclose is used to handle potential floating point inaccuracies.
        dot_product = np.dot(nv_axis, e_field_direction)

        if np.isclose(dot_product, 0):
            # Perpendicular (transverse) case: E is perpendicular to NV axis.
            # This causes splitting into two resonance peaks.
            perpendicular_orientations_count += 1
            resonance_frequencies.add("D + delta")
            resonance_frequencies.add("D - delta")
        else:
            # Parallel (longitudinal) case: E is parallel to NV axis.
            # This does not cause splitting, resulting in one peak.
            parallel_orientations_count += 1
            # Note: The position D_prime might be shifted but it is a single peak.
            resonance_frequencies.add("D_prime")

    # The number of resonances from each group
    resonances_from_parallel_group = 1
    resonances_from_perpendicular_groups = 2 # The two perpendicular groups are equivalent

    # 4. Print the results and the final equation.
    print(f"There are {parallel_orientations_count} orientation group(s) parallel to the E-field.")
    print(f"This group contributes {resonances_from_parallel_group} unique resonance(s).\n")
    
    print(f"There are {perpendicular_orientations_count} orientation group(s) perpendicular to the E-field.")
    print(f"These groups are equivalent and together contribute {resonances_from_perpendicular_groups} unique resonance(s).\n")

    total_resonances = len(resonance_frequencies)

    print("The final equation for the total number of unique resonances is:")
    # We output each number in the final equation.
    print(f"{resonances_from_parallel_group} (from parallel group) + {resonances_from_perpendicular_groups} (from perpendicular groups) = {total_resonances}")
    
    print(f"\nTotal number of ODMR resonances: {total_resonances}")

# Execute the function to get the answer.
calculate_odmr_resonances()
<<<3>>>