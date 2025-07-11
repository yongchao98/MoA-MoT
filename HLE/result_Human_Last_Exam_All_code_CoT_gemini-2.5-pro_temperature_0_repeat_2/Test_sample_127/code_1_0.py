import numpy as np

def calculate_odmr_resonances():
    """
    Calculates the number of ODMR resonances for a hypothetical NV center
    in a cubic lattice under an electric field.
    """
    # Step 1: Define the physical setup.
    # The electric field is applied along a cubic lattice edge, which we define as the x-axis.
    e_field_direction = np.array([1, 0, 0])
    
    # In a cubic lattice, the NV center's symmetry axis can align with the three principal <100> axes.
    nv_orientations = [
        np.array([1, 0, 0]),  # Parallel to x-axis
        np.array([0, 1, 0]),  # Parallel to y-axis
        np.array([0, 0, 1])   # Parallel to z-axis
    ]

    print("Analyzing the ODMR resonances for an NV center in a hypothetical cubic lattice.")
    print(f"Electric field is applied along the direction: {e_field_direction}")
    print("Possible NV center orientations are along the <100> axes (x, y, z).\n")

    # Step 2: Use a set to store the unique types of resonances.
    # We use symbolic names because we are counting peaks, not calculating their exact frequencies.
    unique_resonances = set()
    
    num_parallel_resonances = 0
    num_perpendicular_resonances = 0

    # Step 3: Iterate through each possible NV orientation and determine its resonance signature.
    for i, nv_axis in enumerate(nv_orientations):
        # Check if the NV axis is parallel or perpendicular to the E-field
        # using the dot product.
        dot_product = np.dot(e_field_direction, nv_axis)

        if dot_product == 1:  # Parallel case
            print(f"Case {i+1}: NV axis {nv_axis} is PARALLEL to the E-field.")
            print("   - The m_s=+1 and m_s=-1 states are not split.")
            print("   - This orientation contributes 1 resonance peak.")
            unique_resonances.add('parallel_resonance')
            num_parallel_resonances = 1 # There is only one type of parallel resonance
        
        elif dot_product == 0:  # Perpendicular case
            print(f"Case {i+1}: NV axis {nv_axis} is PERPENDICULAR to the E-field.")
            print("   - The m_s=+1 and m_s=-1 degeneracy is lifted, splitting the resonance.")
            print("   - This orientation contributes 2 resonance peaks.")
            unique_resonances.add('perpendicular_resonance_1')
            unique_resonances.add('perpendicular_resonance_2')
            num_perpendicular_resonances = 2 # There are two types of perpendicular resonances
    
    # Step 4: Summarize and print the final result.
    print("\n--- Summary ---")
    print(f"Number of distinct resonance peaks from parallel orientations: {num_parallel_resonances}")
    print(f"Number of distinct resonance peaks from perpendicular orientations: {num_perpendicular_resonances}")
    
    total_resonances = len(unique_resonances)
    
    print("\nThe final equation for the total number of distinct resonances is:")
    print(f"{num_parallel_resonances} (from parallel) + {num_perpendicular_resonances} (from perpendicular) = {total_resonances}")
    
    print(f"\nTherefore, you will see a total of {total_resonances} resonances in ODMR.")
    
    return total_resonances

if __name__ == '__main__':
    final_answer = calculate_odmr_resonances()
    # The final answer is printed within the function, but we can also use it here.
    # print(f"\n<<< {final_answer} >>>")

calculate_odmr_resonances()
<<<3>>>