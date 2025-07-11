import numpy as np

def calculate_odmr_resonances():
    """
    Calculates the number of ODMR resonances for a hypothetical NV center
    in a cubic lattice with an applied electric field.
    """

    # Step 1: Define the possible NV center orientations.
    # In a simple cubic lattice, the NV axis aligns with the principal axes.
    # There are three distinct, orthogonal groups of orientations.
    # We represent them by unit vectors.
    nv_orientations = {
        "along [1,0,0] axis": np.array([1, 0, 0]),
        "along [0,1,0] axis": np.array([0, 1, 0]),
        "along [0,0,1] axis": np.array([0, 0, 1]),
    }

    # Step 2: Define the applied electric field.
    # The field is parallel to one of the cubic lattice edges.
    # We choose the [0,0,1] direction without loss of generality.
    e_field_vector = np.array([0, 0, 1])
    print(f"Applied electric field is along the direction: {e_field_vector.tolist()}\n")

    # The ODMR frequencies are given by f = D' +/- delta, where D' is the
    # center frequency and delta is the splitting.
    # D' = D + d_parallel * E_parallel
    # delta = d_perp * E_perpendicular
    # Let D be the zero-field splitting frequency and E be the E-field magnitude.
    D = "D"
    d_parallel_E = "d_||*E"
    d_perp_E = "d_perp*E"

    # Use a set to store the unique resonance frequency expressions.
    unique_resonances = set()

    print("Analyzing each group of NV centers:")
    print("-----------------------------------")

    # Step 3: Iterate through each orientation and determine the resonance frequencies.
    for name, nv_axis in nv_orientations.items():
        print(f"For NV centers {name}:")
        
        # E_parallel is proportional to the cosine of the angle between E-field and NV axis.
        e_parallel_component = np.dot(e_field_vector, nv_axis)
        
        # E_perpendicular is proportional to the sine of the angle.
        e_perpendicular_component = np.linalg.norm(np.cross(e_field_vector, nv_axis))

        # Calculate the center frequency and splitting
        center_freq_shift = f"({e_parallel_component:.0f} * {d_parallel_E})"
        splitting = f"({e_perpendicular_component:.0f} * {d_perp_E})"

        # Case 1: NV axis is parallel to the E-field. No splitting occurs.
        if np.isclose(e_perpendicular_component, 0):
            resonance = f"{D} + {center_freq_shift}"
            unique_resonances.add(resonance)
            print(f"  - NV axis is PARALLEL to the E-field.")
            print(f"  - This results in 1 resonance line at the frequency:")
            print(f"    f = {D} + {e_parallel_component:.0f} * {d_parallel_E}")
            
        # Case 2: NV axis is perpendicular to the E-field. Splitting occurs.
        else:
            resonance1 = f"{D} + {center_freq_shift} + {splitting}"
            resonance2 = f"{D} + {center_freq_shift} - {splitting}"
            unique_resonances.add(resonance1)
            unique_resonances.add(resonance2)
            print(f"  - NV axis is PERPENDICULAR to the E-field.")
            print(f"  - This results in 2 resonance lines at the frequencies:")
            # Since e_parallel is 0, the center freq shift term is 0.
            print(f"    f1 = {D} + {e_perpendicular_component:.0f} * {d_perp_E}")
            print(f"    f2 = {D} - {e_perpendicular_component:.0f} * {d_perp_E}")
        print("-" * 35)

    # Step 4: Print the final result.
    total_resonances = len(unique_resonances)
    print("\nSummary:")
    print("Combining the signals from all NV centers, we get a total of {} distinct resonance lines.".format(total_resonances))

if __name__ == '__main__':
    calculate_odmr_resonances()