def solve_ising_susceptibility():
    """
    This function outlines the derivation and prints the final expression for the
    magnetic susceptibility of the Ising model on a sparse random graph.
    """

    print("Step 1: Derive the connected correlation function C_l.")
    print("The perturbation propagates along the path from site l to 0. By linearizing the belief propagation equations, we find the correlation function to be:")
    print("C_l = (1 - m_0^2) * A^l")
    print("Here, m_0 is the magnetization at site 0, and A is the message amplification factor.\n")

    print("Step 2: Substitute C_l into the formula for susceptibility chi.")
    print("The susceptibility is given by the sum:")
    print("chi = beta * sum_{l=1 to inf} [c * (c-1)^(l-1) * C_l]")
    print("Substituting our expression for C_l gives a geometric series:")
    print("chi = beta * c * (1 - m_0^2) * A / (1 - (c-1) * A)\n")

    print("Step 3: Express the result using the provided constant N.")
    print("The note defines a constant N = beta * c * (1 - m_0^2) / (c - 1).")
    print("Using N, the susceptibility can be written compactly as:")
    print("chi = N * (c - 1) * A / (1 - (c - 1) * A)\n")

    print("Step 4: Define the amplification factor A.")
    print("The factor A depends on the system parameters and the cavity magnetization m_c:")
    print("A = tanh(beta * J) * (1 - m_c^2) / (1 - (tanh(beta * J) * m_c)^2)")
    print("where J is the coupling constant, c is the connectivity, beta is the inverse temperature, and m_c is the cavity magnetization, which is determined self-consistently.\n")
    
    print("--- Final Answer ---")
    print("The magnetic susceptibility is:")
    print("chi = (N * (c - 1) * A) / (1 - (c - 1) * A)")
    print("With N = (beta * c * (1 - m_0^2)) / (c - 1) and A defined as above.")

solve_ising_susceptibility()