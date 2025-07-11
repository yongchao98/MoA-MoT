def solve_cpn1_mass_ratio():
    """
    Determines and explains the leading-order asymptotic mass ratio
    for excitations in CP(N-1) models.
    """
    print("This script derives the leading-order asymptotic mass ratio between the lightest solitonic excitation and the subsequent higher excitation in the spectrum of CP(N-1) models as N approaches infinity.")
    print("-" * 70)

    print("\nStep 1: Identification of the Physical States")
    print("The problem refers to 'solitonic excitations'. In the context of CP(N-1) models, these are kink states, which are non-perturbative solutions carrying topological charge. The spectrum of these states is distinct from the fundamental particles of the model.")

    print("\nStep 2: Mass Formula for Kink States")
    print("The mass of a kink state, indexed by an integer k = 1, 2, ..., is given by a semi-classical formula:")
    print("  M_k = C(N) * sin(k * pi / N)")
    print("Here, C(N) is a normalization factor proportional to N, and k corresponds to a topological charge of k/N.")

    print("\nStep 3: Masses of the Lightest Solitons")
    print("The 'lightest solitonic excitation' corresponds to the lowest possible value of k, which is k=1.")
    print("  Mass of lightest soliton (M_1) = C(N) * sin(pi / N)")
    print("\nThe 'subsequent higher excitation' corresponds to the next integer value, k=2.")
    print("  Mass of next soliton (M_2) = C(N) * sin(2 * pi / N)")

    print("\nStep 4: Calculation of the Mass Ratio")
    print("The mass ratio R is the mass of the next soliton divided by the mass of the lightest soliton:")
    print("  R = M_2 / M_1 = [C(N) * sin(2 * pi / N)] / [C(N) * sin(pi / N)]")
    print("The factor C(N) cancels out, leaving:")
    print("  R = sin(2 * pi / N) / sin(pi / N)")

    print("\nStep 5: Simplification using Trigonometry")
    print("We apply the double-angle identity sin(2x) = 2*sin(x)*cos(x), setting x = pi/N:")
    print("  R = [2 * sin(pi / N) * cos(pi / N)] / sin(pi / N)")
    print("After canceling the sin(pi / N) term, the ratio simplifies to:")
    print("  R = 2 * cos(pi / N)")

    print("\nStep 6: Asymptotic Limit for Large N")
    print("The problem asks for the ratio as N approaches infinity. In this limit:")
    print("  lim (N->inf) [pi / N] = 0")
    print("Therefore, the cosine term becomes:")
    print("  lim (N->inf) cos(pi / N) = cos(0) = 1")
    print("Substituting this back into the expression for R gives the asymptotic ratio:")
    print("  R_asymptotic = 2 * 1 = 2")

    print("\n--- Final Equation ---")
    print("In the large N limit, the mass of the k-th soliton, M_k, is proportional to k.")
    print("Therefore, the ratio of the masses of the k=2 and k=1 states is:")
    
    mass_k2 = 2
    mass_k1 = 1
    final_ratio = mass_k2 / mass_k1
    
    print(f"Mass Ratio = (Mass of 2nd Soliton) / (Mass of 1st Soliton) = {mass_k2} / {mass_k1} = {int(final_ratio)}")

# Execute the derivation
solve_cpn1_mass_ratio()