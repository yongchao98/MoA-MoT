import numpy as np

def calculate_odmr_resonances():
    """
    Calculates and explains the number of ODMR resonances for a hypothetical
    NV center in a cubic lattice under an electric field parallel to a lattice edge.

    Assumptions:
    1. The NV center's symmetry axis can align along the three principal <100> axes
       of the cubic lattice: [1,0,0], [0,1,0], and [0,0,1].
    2. The total ODMR spectrum is the sum of spectra from these three distinct
       populations of NV centers.
    3. An electric field perpendicular to the NV axis splits the m_s = +/-1
       degeneracy (Stark effect), leading to two resonance lines. An E-field
       parallel to the axis does not cause splitting.
    """

    # Define symbolic constants for the physics parameters to show the equations.
    # D0 is the zero-field splitting frequency.
    # d_perp is the coupling constant for the transverse Stark effect (splitting).
    # d_para is the coupling constant for the longitudinal Stark effect (shifting).
    # E0 is the magnitude of the applied electric field.
    D0_str = "D₀"
    d_perp_str = "d⊥"
    d_para_str = "d∥"
    E0_str = "E₀"

    # Define the electric field vector. It's parallel to a cubic lattice edge,
    # which we define as the x-axis without loss of generality.
    E_vector_direction = np.array([1, 0, 0])

    # Define the possible orientations of the NV center axis.
    nv_orientations = {
        "Group 1 (axis || x)": np.array([1, 0, 0]),
        "Group 2 (axis || y)": np.array([0, 1, 0]),
        "Group 3 (axis || z)": np.array([0, 0, 1]),
    }

    print("Analyzing the ODMR spectrum for a hypothetical NV center in a cubic lattice.")
    print(f"The electric field is applied along the {list(E_vector_direction)} direction.\n")
    print("There are 3 possible orientations for the NV center's axis: along [1,0,0], [0,1,0], and [0,0,1].")
    print("-" * 65)

    resonance_equations = set()

    for group_name, nv_axis in nv_orientations.items():
        print(f"Processing {group_name} with NV axis along {list(nv_axis)}:")

        # The component of the E-field parallel to the NV axis.
        E_parallel_component = np.dot(E_vector_direction, nv_axis)

        # A non-zero perpendicular component causes splitting.
        # This is non-zero if the E-field is not parallel to the NV axis.
        if np.isclose(np.abs(E_parallel_component), 1.0):
            # E-field is parallel to the NV axis. No splitting.
            # Resonance frequency f = D0 + d_para * E0^2
            shift_term = f"{d_para_str}*({E0_str})²"
            eq = f"{D0_str} + {shift_term}"
            print(f"  - The E-field is PARALLEL to the NV axis.")
            print(f"  - This results in NO splitting, only a frequency shift.")
            print(f"  - Observed Resonance Equation: f = {eq}")
            resonance_equations.add(eq)
        else:
            # E-field is perpendicular to the NV axis. Splitting occurs.
            # The longitudinal shift effect is zero since E_parallel is zero.
            # Resonances are f = D0 ± d_perp * E0
            split_term = f"{d_perp_str}*{E0_str}"
            eq1 = f"{D0_str} + {split_term}"
            eq2 = f"{D0_str} - {split_term}"
            print(f"  - The E-field is PERPENDICULAR to the NV axis.")
            print(f"  - This results in a splitting of the resonance line.")
            print(f"  - Observed Resonance Equations: f = {D0_str} ± {split_term}")
            resonance_equations.add(eq1)
            resonance_equations.add(eq2)
        print("-" * 65)

    print("\nSummary of all unique resonance frequency equations:")
    # Sorting for consistent output order
    for i, eq in enumerate(sorted(list(resonance_equations))):
        print(f"  Line {i+1}: f = {eq}")

    num_resonances = len(resonance_equations)
    print(f"\nIn total, we observe {num_resonances} distinct resonance lines in the ODMR spectrum.")
    return num_resonances

if __name__ == '__main__':
    final_answer = calculate_odmr_resonances()
    print(f"\n<<<{final_answer}>>>")