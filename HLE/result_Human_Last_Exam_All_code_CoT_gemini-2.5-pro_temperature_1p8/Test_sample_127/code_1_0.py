import numpy as np

def calculate_odmr_resonances():
    """
    Calculates the number of ODMR resonances for a hypothetical NV center
    in a cubic lattice under an electric field.
    """
    # 1. Define the physical parameters
    # Zero-field splitting (D) in GHz. This is the central resonance frequency without any splitting.
    D = 2.87
    # Magnitude of the Stark splitting (δ) in GHz, induced by the electric field.
    # We choose an arbitrary non-zero value for demonstration.
    stark_splitting_delta = 0.05

    # 2. Define the system geometry
    # In a cubic lattice, the three high-symmetry orientations are along the cartesian axes.
    nv_orientations = [
        np.array([1, 0, 0]),  # Orientation along X-axis
        np.array([0, 1, 0]),  # Orientation along Y-axis
        np.array([0, 0, 1]),  # Orientation along Z-axis
    ]
    # The electric field is applied along one edge. Let's choose the Z-axis.
    electric_field_direction = np.array([0, 0, 1])

    all_resonances = []
    print("Analyzing each possible NV center orientation:")
    print("-" * 50)

    # 3. Iterate through each orientation and calculate resonances
    for i, nv_axis in enumerate(nv_orientations):
        # Calculate the cosine of the angle between E-field and NV axis
        cos_theta = np.dot(nv_axis, electric_field_direction)

        # Check if the vectors are parallel or perpendicular
        # np.isclose is used to handle potential floating point inaccuracies.
        if np.isclose(np.abs(cos_theta), 1.0):
            # E-field is PARALLEL to the NV axis.
            # Perpendicular component of E is zero, so no Stark splitting.
            print(f"Orientation {i+1} ({str(list(nv_axis))}): Parallel to E-field.")
            print(f"  - Result: No splitting. One resonance at D.")
            print(f"  - Frequency = {D:.2f} GHz\n")
            all_resonances.append(D)
        elif np.isclose(cos_theta, 0.0):
            # E-field is PERPENDICULAR to the NV axis.
            # Perpendicular component of E is maximal, causing splitting.
            freq_minus = D - stark_splitting_delta
            freq_plus = D + stark_splitting_delta
            print(f"Orientation {i+1} ({str(list(nv_axis))}): Perpendicular to E-field.")
            print(f"  - Result: States m_s=±1 split. Two resonances at D±δ.")
            print(f"  - Frequencies = {D:.2f} - {stark_splitting_delta:.2f} = {freq_minus:.2f} GHz")
            print(f"                  {D:.2f} + {stark_splitting_delta:.2f} = {freq_plus:.2f} GHz\n")
            all_resonances.append(freq_minus)
            all_resonances.append(freq_plus)

    # 4. Count the number of unique resonances
    # Using a set automatically handles duplicates.
    unique_resonances = sorted(list(set(all_resonances)))
    
    print("-" * 50)
    print("Summary:")
    print(f"The combined ODMR spectrum will show resonances at the following unique frequencies:")
    print(f"{[f'{freq:.2f} GHz' for freq in unique_resonances]}")
    print(f"\nThe total number of distinct resonances is: {len(unique_resonances)}")


if __name__ == "__main__":
    calculate_odmr_resonances()