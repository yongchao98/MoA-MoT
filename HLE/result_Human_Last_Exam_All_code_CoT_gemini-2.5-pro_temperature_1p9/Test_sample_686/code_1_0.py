def solve_ising_susceptibility():
    """
    This function provides a step-by-step derivation of the magnetic susceptibility
    for an Ising model on a sparse random graph and prints the final result.
    """

    print("### Derivation of Magnetic Susceptibility ###\n")

    # Step 1: State the initial formulas
    print("Step 1: Initial formulas")
    print("--------------------------")
    print("The magnetic susceptibility is given by:")
    print("  chi = beta * sum_{l=1 to inf} [c * (c-1)^(l-1) * C_l]  (1)")
    print("\nThe connected correlation C_l is related to the derivative of the magnetization m_0:")
    print("  C_l = (1/beta) * dm_0 / dB_l  (2)")
    print("\nWhere m_0 = <sigma_0> is the magnetization at site 0, and B_l is the field at site l.")
    print("--------------------------\n")

    # Step 2: Model the propagation of the field perturbation
    print("Step 2: Calculate the correlation function C_l")
    print("--------------------------")
    print("On a sparse random graph (locally a tree), a perturbation delta_B_l at site l propagates")
    print("along the unique path to site 0. Let's trace this effect.")
    print("The magnetization at site 0 is m_0 = tanh[beta*(B_0 + sum_{j in neighbors} u_{j->0})].")
    print("The derivative with respect to B_l is:")
    print("  dm_0 / dB_l = beta * (1 - m_0^2) * (d u_{1->0} / dB_l)")
    print("where site 1 is the neighbor of 0 on the path from l.")
    print("\nThe message propagation can be written as a recursion:")
    print("  d u_{i-1 -> i-2} / dB_l = K * (d u_{i -> i-1} / dB_l)")
    print("where K is the propagation factor.")
    print("\nThis recursion unwraps to:")
    print("  d u_{1->0} / dB_l = K^(l-1) * (d u_{l->l-1} / dB_l)")
    print("\nThe initial perturbation step is found to be:")
    print("  d u_{l->l-1} / dB_l = K")
    print("\nTherefore, d u_{1->0} / dB_l = K^l.")
    print("\nSubstituting back, we get dm_0 / dB_l = beta * (1 - m_0^2) * K^l.")
    print("\nUsing equation (2), we find the correlation function:")
    print("  C_l = (1/beta) * [beta * (1 - m_0^2) * K^l] = (1 - m_0^2) * K^l  (3)")
    print("--------------------------\n")

    # Step 3: Substitute C_l into the sum for chi
    print("Step 3: Evaluate the sum for susceptibility chi")
    print("--------------------------")
    print("Substituting C_l from (3) into equation (1):")
    print("  chi = beta * sum_{l=1 to inf} [c * (c-1)^(l-1) * (1 - m_0^2) * K^l]")
    print("\nRearranging the terms:")
    print("  chi = beta * c * (1 - m_0^2) * sum_{l=1 to inf} [(c-1)^(l-1) * K^l]")
    print("  chi = beta * [c / (c-1)] * (1 - m_0^2) * sum_{l=1 to inf} [((c-1) * K)^l]")
    print("\nThe sum is a geometric series sum_{l=1 to inf} r^l = r / (1 - r) with r = (c-1)*K.")
    print("This requires |(c-1)*K| < 1 for convergence.")
    print("\nEvaluating the sum gives:")
    print("  sum = [(c-1)*K] / [1 - (c-1)*K]")
    print("\nSo the susceptibility is:")
    print("  chi = beta * [c / (c-1)] * (1 - m_0^2) * {[(c-1)*K] / [1 - (c-1)*K]}")
    print("  chi = beta * c * (1 - m_0^2) * K / [1 - (c-1)*K]")
    print("--------------------------\n")
    
    # Step 4: Final expression
    print("Step 4: Final expression for chi")
    print("--------------------------")
    print("Using the provided abbreviation N = beta * c * (1 - m_0^2) / (c-1), we can rewrite chi.")
    print("  chi = N * (c-1) * K / [1 - (c-1) * K]")
    print("\nThis is the final expression for the magnetic susceptibility.")
    print("\nFinal Equation Breakdown:")
    print("chi = N * (c - 1) * K / (1 - (c - 1) * K)")
    print("Where:")
    print("  chi: Magnetic susceptibility")
    print("  N: A constant defined as beta * c * (1 - m_0^2) / (c - 1)")
    print("  beta: Inverse temperature (1/(k_B*T))")
    print("  c: Connectivity of the random graph (number of neighbors per site)")
    print("  m_0: Magnetization per site in the homogeneous state")
    print("  K: The propagation factor, defined as K = tanh(beta*J) * (1 - m_c^2) / (1 - (tanh(beta*J)*m_c)^2)")
    print("  J: The coupling constant between spins")
    print("  m_c: The cavity magnetization, found by solving the self-consistency equation m_c = tanh(beta*B + (c-1)*arctanh(tanh(beta*J)*m_c))")

solve_ising_susceptibility()

print("\n<<<chi = N * (c-1) * K / (1 - (c-1)*K)>>>")