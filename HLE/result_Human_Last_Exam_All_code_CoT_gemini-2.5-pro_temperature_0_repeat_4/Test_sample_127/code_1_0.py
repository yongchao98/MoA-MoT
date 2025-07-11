import numpy as np

def calculate_odmr_resonances():
    """
    Calculates the number of ODMR resonances for a hypothetical NV center
    in a cubic lattice under an electric field.
    """
    # Step 1: Define the possible orientations for the NV center in a cubic lattice.
    # These are along the principal axes <100>.
    nv_orientations = [
        np.array([1, 0, 0]),  # Along x-axis
        np.array([0, 1, 0]),  # Along y-axis
        np.array([0, 0, 1]),  # Along z-axis
    ]

    # Step 2: Define the applied electric field, parallel to one edge (e.g., x-axis).
    e_field_direction = np.array([1, 0, 0])

    print("System Configuration:")
    print(f"  - Hypothetical NV center in a CUBIC lattice.")
    print(f"  - Possible NV orientations: {[1,0,0]}, {[0,1,0]}, {[0,0,1]}")
    print(f"  - Applied Electric Field direction: {e_field_direction.tolist()}")
    print("-" * 50)

    # Step 3: Analyze the effect of the E-field on each orientation.
    num_parallel_orientations = 0
    num_perpendicular_orientations = 0

    print("Analysis of each NV orientation:")
    for i, nv_axis in enumerate(nv_orientations):
        # Use the dot product to check for parallelism/perpendicularity.
        # For unit vectors, a dot product of 1 means parallel, 0 means perpendicular.
        dot_product = np.dot(nv_axis, e_field_direction)

        if np.isclose(np.abs(dot_product), 1):
            print(f"  - Orientation {nv_axis.tolist()}: PARALLEL to E-field. Contributes 1 resonance line.")
            num_parallel_orientations += 1
        elif np.isclose(dot_product, 0):
            print(f"  - Orientation {nv_axis.tolist()}: PERPENDICULAR to E-field. Contributes 2 resonance lines.")
            num_perpendicular_orientations += 1

    print("-" * 50)

    # Step 4: Calculate the total number of unique resonance lines.
    # The single parallel orientation gives one unique line.
    # The two perpendicular orientations give the same pair of two unique lines.
    resonances_from_parallel = 1 if num_parallel_orientations > 0 else 0
    resonances_from_perpendicular = 2 if num_perpendicular_orientations > 0 else 0

    total_resonances = resonances_from_parallel + resonances_from_perpendicular

    print("Final Calculation:")
    print(f"Number of resonance lines from the parallel group = {resonances_from_parallel}")
    print(f"Number of resonance lines from the perpendicular group = {resonances_from_perpendicular}")
    print("\nThe total number of distinct ODMR resonances is the sum of unique lines from each group.")
    print(f"Total Resonances = {resonances_from_parallel} + {resonances_from_perpendicular} = {total_resonances}")

calculate_odmr_resonances()
<<<3>>>