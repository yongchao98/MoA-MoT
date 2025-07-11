def solve_mass_ratio():
    """
    Calculates the mass ratio for solitonic excitations in the CP(N-1) model.

    In the large-N limit, the mass M_k of a BPS solitonic state with
    topological charge k is directly proportional to k.
    M_k = k * M_1, where M_1 is the mass of the fundamental soliton (k=1).
    """

    # The lightest solitonic excitation corresponds to a topological charge of k=1.
    k_lightest = 1

    # The subsequent higher excitation corresponds to a state with mass
    # equivalent to two fundamental solitons, which can be a bound state with k=2.
    k_next = 2

    # We are asked for the ratio of the mass of the subsequent excitation (M_2)
    # to the mass of the lightest solitonic excitation (M_1).
    # Ratio = M_2 / M_1 = (k_next * M_1) / (k_lightest * M_1)
    
    ratio = k_next / k_lightest

    print("This script calculates the leading-order asymptotic mass ratio between the lightest solitonic excitation and the subsequent higher excitation in the spectrum of CP(N-1) models.")
    print("-" * 50)
    print("The mass of a solitonic state with topological charge k is M_k.")
    print("The relationship is M_k = k * M_1, where M_1 is the mass of the lightest soliton (k=1).\n")

    print("1. The lightest solitonic excitation has k=1, so its mass is M_1.")
    print(f"   - Let k_lightest = {k_lightest}")

    print("\n2. The subsequent higher excitation has a mass equivalent to two lightest solitons, so its mass is M_2.")
    print(f"   - Let k_next = {k_next}")
    print("   - M_2 = 2 * M_1\n")

    print("3. The mass ratio is M_2 / M_1.")
    print("\nFinal Equation:")
    print(f"Ratio = M_{k_next} / M_{k_lightest} = ({k_next} * M_1) / ({k_lightest} * M_1) = {int(ratio)}")

solve_mass_ratio()
<<<2>>>