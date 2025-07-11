def solve_susceptibility():
    """
    This function prints the derived formula for the magnetic susceptibility chi.
    The derivation is based on the cavity method for an Ising model on a sparse random graph.
    """
    
    # The derived magnetic susceptibility is given by:
    # chi = N * (c-1) * T_u / (1 - (c-1)*T_u)
    # where N is a given constant, c is the connectivity,
    # and T_u is the propagation factor for a perturbation on the graph.
    # T_u depends on the temperature, coupling strength, and the spontaneous magnetization.

    print("The magnetic susceptibility chi is found to be:")
    print("\u03C7 = N * ( (c - 1) * T_u ) / ( 1 - (c - 1) * T_u )")
    print("\nwhere:")
    print("N = \u03B2 * c * (1 - m\u2080\u00B2) / (c - 1) is a given constant.")
    print("c is the connectivity of the graph.")
    print("m\u2080 is the spontaneous magnetization per site.")
    print("\u03B2 is the inverse temperature.")
    print("T_u is the perturbation propagation factor on the Bethe lattice.")

solve_susceptibility()
