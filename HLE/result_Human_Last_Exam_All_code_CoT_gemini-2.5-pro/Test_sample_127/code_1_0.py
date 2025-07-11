import numpy as np

def calculate_odmr_resonances():
    """
    Calculates the number of ODMR resonances for a hypothetical NV center
    in a cubic lattice under an electric field.
    """
    # Define constants for the simulation. The exact values do not change the
    # number of lines, but they help illustrate the calculation.
    D_ZFS = 2.870  # Zero-field splitting in GHz
    E_FIELD_MAGNITUDE = 1.0  # Arbitrary units for electric field strength
    STARK_CONSTANT = 0.1   # Stark shift constant in GHz / (E-field unit)

    # In a cubic lattice, there are three equivalent principal axes.
    # These are the possible orientations for our hypothetical NV center.
    orientations = {
        "Axis [1,0,0]": np.array([1, 0, 0]),
        "Axis [0,1,0]": np.array([0, 1, 0]),
        "Axis [0,0,1]": np.array([0, 0, 1]),
    }

    # The electric field is applied along one of the cubic lattice edges.
    # We choose the [100] direction without loss of generality.
    e_field_vector = np.array([E_FIELD_MAGNITUDE, 0, 0])

    print("--- ODMR Resonance Calculation ---")
    print(f"Scenario: NV-like center in a cubic lattice.")
    print(f"Zero-Field Splitting (D) = {D_ZFS:.3f} GHz")
    print(f"Electric Field (E) applied along [1,0,0]\n")

    all_resonances = set()

    # Analyze the effect of the E-field on each NV orientation
    for name, nv_axis in orientations.items():
        print(f"Analyzing Population with {name}:")

        # The Stark shift is proportional to the component of E perpendicular to the NV axis.
        # E_perp = E - (E . nv_axis) * nv_axis
        e_dot_nv = np.dot(e_field_vector, nv_axis)
        e_parallel = e_dot_nv * nv_axis
        e_perpendicular_vector = e_field_vector - e_parallel
        e_perp_magnitude = np.linalg.norm(e_perpendicular_vector)

        # Calculate the energy shift (delta)
        stark_shift = STARK_CONSTANT * e_perp_magnitude

        # If the shift is (close to) zero, there is no splitting
        if np.isclose(stark_shift, 0):
            print(f"  - E-field is parallel to the NV axis.")
            print(f"  - Stark shift is 0.")
            print(f"  - A single resonance is observed at D.")
            print(f"  - Equation: {D_ZFS:.3f} GHz")
            all_resonances.add(round(D_ZFS, 3))
        # If the shift is non-zero, the m_s=+/-1 states split
        else:
            res1 = D_ZFS - stark_shift
            res2 = D_ZFS + stark_shift
            print(f"  - E-field is perpendicular to the NV axis.")
            print(f"  - Stark shift (δ) is {stark_shift:.3f} GHz.")
            print(f"  - Two resonances are observed at D ± δ.")
            print(f"  - Equations: {D_ZFS:.3f} - {stark_shift:.3f} = {res1:.3f} GHz")
            print(f"                 {D_ZFS:.3f} + {stark_shift:.3f} = {res2:.3f} GHz")
            all_resonances.add(round(res1, 3))
            all_resonances.add(round(res2, 3))
        print("-" * 35)

    print("\n--- Summary ---")
    print(f"The set of unique resonance frequencies is: {sorted(list(all_resonances))} GHz")
    print(f"Total number of distinct resonances observed = {len(all_resonances)}")

if __name__ == '__main__':
    calculate_odmr_resonances()
<<<3>>>