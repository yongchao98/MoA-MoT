def calculate_mass_ratio():
    """
    Calculates the leading-order asymptotic mass ratio between the lightest
    and subsequent higher solitonic excitations in the CP(N-1) model.

    In the large-N limit, the mass of a solitonic excitation (baryon) with
    topological charge k is given by M_k = |k| * M_1, where M_1 is the mass
    of the fundamental k=1 soliton.
    """

    # The lightest solitonic excitation has topological charge k=1.
    k_lightest = 1

    # The subsequent higher solitonic excitation has topological charge k=2.
    k_subsequent = 2

    # The mass ratio M_2 / M_1 is therefore equal to the ratio of their charges.
    mass_ratio = float(k_subsequent) / float(k_lightest)

    print("Determining the mass ratio of solitonic excitations in the CP(N-1) model.")
    print("Let M_k be the mass of a soliton with topological charge k.")
    print(f"The lightest solitonic excitation has charge k = {k_lightest}, with mass M_{k_lightest}.")
    print(f"The subsequent higher solitonic excitation has charge k = {k_subsequent}, with mass M_{k_subsequent}.")
    print("\nAt leading order in the large-N expansion, the mass ratio is the ratio of the topological charges:")
    print(f"Mass Ratio = M_{k_subsequent} / M_{k_lightest} = {k_subsequent} / {k_lightest} = {mass_ratio}")

calculate_mass_ratio()
<<<2>>>