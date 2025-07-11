import numpy as np

def solve_hypothetical_nv_odmr():
    """
    Calculates the number of ODMR resonances for a hypothetical NV center in a
    cubic lattice under an electric field parallel to a lattice edge.
    """
    
    # --- Step 1 & 2: Define the physical model and E-field ---
    print("This script calculates the number of ODMR resonances for an NV center in a hypothetical cubic lattice.")
    print("Here is the step-by-step analysis:\n")
    print("1. In a cubic lattice, there are three equivalent high-symmetry axes <100>.")
    print("   We assume the NV centers align along these three directions: [1,0,0], [0,1,0], and [0,0,1].")
    print("2. An electric field 'E' is applied along one edge. We'll align it with the z-axis, [0,0,1].")
    print("3. The resonance frequency shift is Δf = (d_par * E_parallel) +/- (d_perp * E_perp).")
    print("-" * 70)

    # The three possible NV center orientation vectors
    nv_orientations = {
        "NVs along [1,0,0] (x-axis)": np.array([1, 0, 0]),
        "NVs along [0,1,0] (y-axis)": np.array([0, 1, 0]),
        "NVs along [0,0,1] (z-axis)": np.array([0, 0, 1])
    }

    # The electric field is along the z-axis. We use a unit vector for direction.
    # The magnitude is represented by 'E' in the final equations.
    e_field_direction = np.array([0, 0, 1])

    # --- Step 3: Calculate the Stark effect for each orientation ---
    # We will use a set to store the unique resonance equations.
    # Each equation will be represented by a tuple of coefficients (coeff_par, coeff_perp).
    unique_resonance_coeffs = set()

    print("Analyzing each group of NV centers:\n")
    for name, nv_axis in nv_orientations.items():
        print(f"For {name}:")
        
        # Calculate the projection of E onto the NV axis to find the parallel component coefficient.
        coeff_par = np.dot(e_field_direction, nv_axis)
        
        # Calculate the perpendicular component vector and its magnitude coefficient.
        e_perp_vector = e_field_direction - coeff_par * nv_axis
        coeff_perp = np.linalg.norm(e_perp_vector)

        print(f"  - E-field parallel to NV axis: {coeff_par:.1f} * E")
        print(f"  - E-field perpendicular to NV axis: {coeff_perp:.1f} * E")

        if np.isclose(coeff_perp, 0):
            # No perpendicular component, so no splitting. One resonance.
            print("  - Result: No splitting. One resonance line emerges.\n")
            unique_resonance_coeffs.add((coeff_par, 0.0))
        else:
            # Perpendicular component causes splitting. Two resonances.
            print("  - Result: Splitting occurs. Two resonance lines emerge.\n")
            unique_resonance_coeffs.add((coeff_par, coeff_perp))
            unique_resonance_coeffs.add((coeff_par, -coeff_perp))

    # --- Step 4: Count and display the final unique resonances ---
    print("-" * 70)
    print(f"By combining the results from all NV orientations, we find a total of {len(unique_resonance_coeffs)} unique resonance lines.")
    print("\nThe final equations for the resonance frequencies are (where 'D' is the zero-field splitting):")
    
    # Sort for a consistent and logical output order
    sorted_coeffs = sorted(list(unique_resonance_coeffs), key=lambda c: (c[0], -c[1]))
    
    for i, (c_par, c_perp) in enumerate(sorted_coeffs):
        parts = ["f = D"]
        if not np.isclose(c_par, 0):
            parts.append(f"{c_par:+.1f}*d_par*E")
        if not np.isclose(c_perp, 0):
            parts.append(f"{c_perp:+.1f}*d_perp*E")
        
        # Join the parts into a clean string, e.g., "f = D + 1.0*d_par*E"
        equation_str = " ".join(parts).replace(" +", " + ").replace(" -", " - ")
        print(f"  {i+1}. {equation_str}")

# Execute the function to get the answer
solve_hypothetical_nv_odmr()