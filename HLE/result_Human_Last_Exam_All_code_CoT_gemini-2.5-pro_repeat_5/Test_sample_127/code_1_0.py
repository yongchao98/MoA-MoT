import numpy as np

def calculate_odmr_resonances():
    """
    Calculates the number of ODMR resonances for a hypothetical NV center
    in a cubic lattice under an electric field parallel to a lattice edge.
    """

    # 1. Define the possible NV orientations in a cubic lattice.
    # These are the three principal, orthogonal directions.
    nv_orientations = {
        "orientation_1": np.array([1, 0, 0]),
        "orientation_2": np.array([0, 1, 0]),
        "orientation_3": np.array([0, 0, 1])
    }

    # 2. Define the direction of the applied electric field.
    # We align it with one of the lattice edges without loss of generality.
    electric_field_direction = np.array([1, 0, 0])

    print("Analyzing the ODMR spectrum for a hypothetical NV center in a cubic lattice.")
    print(f"Applied electric field is along the {electric_field_direction} direction.\n")

    resonance_energies = set()
    num_parallel_resonances = 0
    num_perpendicular_resonances = 0

    # 3. Analyze the effect on each NV orientation.
    # The Stark effect splits the m_s=+/-1 degeneracy if E is perpendicular to the NV axis.
    # Let D be the zero-field splitting and d_perp be the transverse Stark shift parameter.
    # The splitting magnitude is proportional to the perpendicular E-field component.
    
    # Case 1: NV axis is parallel to the E-field.
    # This corresponds to orientation_1.
    nv_axis_parallel = nv_orientations["orientation_1"]
    # The component of E perpendicular to the NV axis is zero.
    # cos(theta) = 1, sin(theta) = 0.
    # No splitting of the m_s=+/-1 states occurs.
    num_parallel_resonances = 1
    print(f"For NV centers with their axis parallel to the E-field ({nv_axis_parallel}):")
    print("The perpendicular E-field component is zero. No splitting occurs.")
    print(f"Result: {num_parallel_resonances} resonance line.\n")

    # Case 2: NV axis is perpendicular to the E-field.
    # This corresponds to orientation_2 and orientation_3.
    # The component of E perpendicular to the NV axis is non-zero.
    # cos(theta) = 0, sin(theta) = 1.
    # The m_s=+/-1 states are split into two new states.
    num_perpendicular_resonances = 2
    nv_axis_perp_1 = nv_orientations["orientation_2"]
    nv_axis_perp_2 = nv_orientations["orientation_3"]
    print(f"For NV centers with their axis perpendicular to the E-field ({nv_axis_perp_1} and {nv_axis_perp_2}):")
    print("The perpendicular E-field component is non-zero, causing a Stark shift.")
    print("The single resonance splits into two distinct lines.")
    print(f"Result: {num_perpendicular_resonances} resonance lines.\n")
    
    # 4. Sum the unique resonance lines.
    # The lines from the perpendicular group are distinct from the line from the parallel group.
    # The two perpendicular orientations are equivalent and their resonances will overlap perfectly.
    total_resonances = num_parallel_resonances + num_perpendicular_resonances
    
    print("--- Summary ---")
    print("We have one group of NV centers parallel to the E-field and one group perpendicular to it.")
    print(f"The final equation for the total number of distinct resonances is:")
    print(f"{num_parallel_resonances} (from parallel NVs) + {num_perpendicular_resonances} (from perpendicular NVs) = {total_resonances}")

calculate_odmr_resonances()