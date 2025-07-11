import numpy as np

def calculate_odmr_resonances():
    """
    Calculates the number of ODMR resonances for a hypothetical NV center in a
    cubic lattice under an electric field parallel to a lattice edge.
    """
    # 1. Define the possible NV orientations in a cubic lattice.
    # These align with the primary axes [100], [010], [001].
    orientations = [
        np.array([1, 0, 0]),
        np.array([0, 1, 0]),
        np.array([0, 0, 1])
    ]

    # 2. Define the electric field vector.
    # It's applied along one of the edges, let's say the z-axis.
    # We use a unit vector as only the direction matters for this analysis.
    e_field = np.array([0, 0, 1])

    print("Analyzing the Stark effect on different NV orientations:")
    print(f"Electric field direction: {e_field}")
    print("-" * 50)

    # 3. Calculate the splitting for each orientation.
    # The splitting is proportional to the magnitude of the E-field component
    # perpendicular to the NV axis. We can use the squared magnitude to avoid
    # square roots, as we only need to distinguish between zero and non-zero splitting.
    
    squared_splitting_magnitudes = []
    for i, n_vector in enumerate(orientations):
        # Project E-field onto the NV axis to find the parallel component
        e_parallel_component = np.dot(e_field, n_vector) * n_vector
        
        # Subtract the parallel component from the total E-field to get the perpendicular component
        e_perpendicular_component = e_field - e_parallel_component
        
        # The magnitude of the splitting is proportional to |E_perp|
        # We calculate the squared magnitude for simplicity.
        squared_mag = np.sum(e_perpendicular_component**2)
        squared_splitting_magnitudes.append(round(squared_mag, 5)) # Round to handle potential float precision issues
        
        print(f"NV Orientation {n_vector}:")
        print(f"  - E-field component parallel to NV axis: {e_parallel_component}")
        print(f"  - E-field component perpendicular to NV axis: {e_perpendicular_component}")
        if squared_mag == 0:
            print("  - Result: No splitting (1 resonance peak)")
        else:
            print("  - Result: Splitting into 2 resonance peaks")
        print()

    # 4. Count the number of unique resonance groups.
    # A set will store only the unique values of splitting magnitudes.
    unique_splittings = set(squared_splitting_magnitudes)
    
    num_unsplit_groups = 0
    if 0.0 in unique_splittings:
        num_unsplit_groups = 1
        
    num_split_groups = len(unique_splittings) - num_unsplit_groups
    
    # Each unsplit group contributes 1 peak.
    # Each split group contributes 2 peaks.
    total_resonances = (num_unsplit_groups * 1) + (num_split_groups * 2)

    print("-" * 50)
    print("Final Calculation:")
    print(f"There is {num_unsplit_groups} group of NV centers with zero splitting.")
    print(f"There are {num_split_groups} distinct groups of NV centers with non-zero splitting.")
    print(f"\nTotal number of resonances = ({num_unsplit_groups} * 1) + ({num_split_groups} * 2)")
    print(f"                             = {num_unsplit_groups * 1} + {num_split_groups * 2}")
    print(f"                             = {total_resonances}")
    print("-" * 50)


calculate_odmr_resonances()

# The final answer is the total number of resonances calculated.
# The code above prints this result along with the step-by-step logic.
# Based on the logic, the number of resonances is 3.
print("\n<<<3>>>")