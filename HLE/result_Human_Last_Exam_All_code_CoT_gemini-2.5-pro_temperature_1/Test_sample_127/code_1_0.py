import numpy as np

def calculate_cubic_nv_odmr():
    """
    Calculates the number of ODMR resonances for a hypothetical NV center
    in a cubic lattice under an electric field parallel to a lattice edge.
    """
    # 1. Define the possible orientations for an NV-like center in a cubic lattice.
    # These are along the principal axes of the cube.
    nv_orientations = {
        "x-axis": np.array([1, 0, 0]),
        "y-axis": np.array([0, 1, 0]),
        "z-axis": np.array([0, 0, 1]),
    }

    # 2. Define the applied electric field. The problem states it is parallel to a
    # cubic lattice edge. We can choose the z-axis without loss of generality.
    # We use a magnitude of 1.0 in arbitrary units for the calculation.
    E_field_vec = np.array([0, 0, 1.0])
    E_mag = np.linalg.norm(E_field_vec)

    # 3. Define placeholder physical constants for the calculation.
    # D is the zero-field splitting frequency.
    # d_perp is the Stark coefficient for the perpendicular electric field.
    # The splitting of the m_s=+-1 states is 2*delta, where delta = d_perp * |E_perp|.
    D_val = 2.87  # in GHz
    d_perp_val = 0.05 # in GHz / (arb. unit of E)

    print("Analyzing a hypothetical NV center in a cubic lattice.")
    print(f"An electric field is applied along the {E_field_vec} direction.")
    print("-" * 50)

    # Use a set to store the unique resonance frequencies found.
    # We store them as tuples (symbolic_string, numeric_value) to handle degeneracy.
    unique_resonances = set()

    # 4. Iterate through each possible NV orientation and calculate its resonance(s).
    for name, nv_axis_vec in nv_orientations.items():
        print(f"Processing NV centers oriented along the {name} ({nv_axis_vec})...")

        # The magnitude of the E-field component perpendicular to the NV axis
        # determines the splitting.
        # E_parallel_component = E_field_vec (dot) nv_axis_vec
        # E_perp_mag = sqrt(|E|^2 - |E_parallel_component|^2)
        E_parallel_comp = np.dot(E_field_vec, nv_axis_vec)
        E_perp_mag = np.sqrt(E_mag**2 - E_parallel_comp**2)

        # Calculate the Stark shift amount 'delta'.
        delta = d_perp_val * E_perp_mag

        if np.isclose(delta, 0):
            print(f"  - E-field is parallel to the NV axis. Perpendicular component is 0.")
            print(f"  - No Stark splitting occurs.")
            print(f"  - One resonance is observed at the zero-field splitting frequency D.")
            # Equation: f = D
            # Numeric:  f = 2.87
            equation = f"f = D = {D_val:.2f} GHz"
            print(f"    Equation: {equation}\n")
            unique_resonances.add(equation)
        else:
            print(f"  - E-field is perpendicular to the NV axis. Perpendicular component is {E_perp_mag:.1f}*|E|.")
            print(f"  - Stark effect splits the resonance into two lines.")
            # Upper resonance
            freq_high = D_val + delta
            # Equation: f = D + d_perp * |E_perp|
            # Numeric:  f = 2.87 + 0.05
            equation_high = f"f = D + d_perp*|E| = {D_val:.2f} + {delta:.2f} = {freq_high:.2f} GHz"
            print(f"    Equation 1: {equation_high}")
            unique_resonances.add(equation_high)

            # Lower resonance
            freq_low = D_val - delta
            # Equation: f = D - d_perp * |E_perp|
            # Numeric:  f = 2.87 - 0.05
            equation_low = f"f = D - d_perp*|E| = {D_val:.2f} - {delta:.2f} = {freq_low:.2f} GHz"
            print(f"    Equation 2: {equation_low}\n")
            unique_resonances.add(equation_low)

    # 5. Count the number of unique resonances.
    num_resonances = len(unique_resonances)
    print("-" * 50)
    print("Summary:")
    print("The set of unique ODMR resonances observed from all orientations is:")
    # Sort for consistent display order
    for resonance in sorted(list(unique_resonances)):
        print(f"  - {resonance}")
        
    print(f"\nThe total number of distinct resonances is {num_resonances}.")

if __name__ == "__main__":
    calculate_cubic_nv_odmr()