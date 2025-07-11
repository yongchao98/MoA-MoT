import numpy as np

def count_odmr_resonances():
    """
    Calculates the number of ODMR resonances for a hypothetical NV center
    in a cubic lattice under an electric field parallel to a lattice edge.
    """
    # 1. Define the possible orientations of the NV center's axis in a cubic lattice.
    # These are the three <100> directions.
    nv_orientations = [
        np.array([1, 0, 0]),  # along x-axis
        np.array([0, 1, 0]),  # along y-axis
        np.array([0, 0, 1])   # along z-axis
    ]

    # 2. Define the direction of the applied electric field.
    # It's parallel to one edge, so we choose the z-axis without loss of generality.
    e_field_direction = np.array([0, 0, 1])
    
    print("Scenario: A hypothetical NV center in a cubic lattice.")
    print(f"Electric field is applied along the {e_field_direction} direction.\n")
    print("Analyzing the effect on each possible NV orientation:")
    print("-" * 60)

    num_parallel_groups = 0
    num_perpendicular_groups = 0

    # 3. Analyze the interaction for each orientation.
    for i, orientation in enumerate(nv_orientations):
        # Use the dot product to determine if the vectors are parallel or perpendicular.
        dot_product = np.dot(orientation, e_field_direction)

        if dot_product == 1:
            # E-field is PARALLEL to the NV axis
            relationship = "parallel"
            resonances = 1
            num_parallel_groups += 1
            # Store the number of resonances for the final equation
            resonances_from_parallel_case = resonances
        elif dot_product == 0:
            # E-field is PERPENDICULAR to the NV axis
            relationship = "perpendicular"
            resonances = 2
            num_perpendicular_groups += 1
            # Store the number of resonances for the final equation
            # All perpendicular groups are equivalent and produce the same doublet.
            resonances_from_perpendicular_case = resonances
        else:
            # This case shouldn't be reached with <100> axes
            relationship = "at an angle"
            resonances = "unknown"

        print(f"NV Group {i+1} (axis along {list(orientation)}):")
        print(f"  - Relationship to E-field: {relationship}")
        print(f"  - Result: The degeneracy of m_s=±1 is NOT lifted." if relationship == "parallel" else "  - Result: The degeneracy of m_s=±1 IS lifted (Stark effect).")
        print(f"  - Contribution: {resonances} resonance line(s).\n")

    # 4. Count the total number of distinct resonance lines.
    # The two perpendicular groups produce the same set of 2 lines.
    # The one parallel group produces 1 unique line.
    total_distinct_resonances = resonances_from_parallel_case + resonances_from_perpendicular_case

    print("-" * 60)
    print("Summary:")
    print(f"The {num_perpendicular_groups} groups of NV centers with axes perpendicular to the E-field are equivalent.")
    print("They collectively produce one set of split resonance lines.\n")

    print("The final calculation for the total number of distinct resonance lines is:")
    print(f"Number of lines from perpendicular case + Number of lines from parallel case")
    print(f"{resonances_from_perpendicular_case} + {resonances_from_parallel_case} = {total_distinct_resonances}")
    
# Execute the function
if __name__ == "__main__":
    count_odmr_resonances()
    
    # Final answer based on the logic
    # print("\n<<<3>>>") # The final answer is enclosed in the required format
