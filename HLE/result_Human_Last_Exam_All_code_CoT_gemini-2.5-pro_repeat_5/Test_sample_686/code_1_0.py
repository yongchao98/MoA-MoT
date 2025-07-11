def print_susceptibility_formula():
    """
    Prints the derived formula for the magnetic susceptibility chi.
    """
    
    # The final expression for the susceptibility chi is composed of several parts.
    
    # 1. The main formula for chi in terms of the constant N, connectivity c, and propagation factor lambda.
    print("The magnetic susceptibility chi is given by the formula:")
    print("chi = N * ( (c - 1) * lambda ) / ( 1 - (c - 1) * lambda )")
    print("-" * 20)
    
    # 2. The definition of the constant N, as provided in the problem.
    print("The constant N is defined as:")
    print("N = beta * c * (1 - m_0**2) / (c - 1)")
    print("-" * 20)

    # 3. The definition of the propagation factor lambda.
    print("The factor 'lambda' represents how a perturbation propagates along the graph. It is defined as:")
    print("lambda = (tanh(beta * J)**2 - tanh(beta * u)**2) / (tanh(beta * J) * (1 - tanh(beta * u)**2))")
    print("-" * 20)

    # 4. The self-consistency equations that determine the values of the message 'u' and magnetization 'm_0'
    #    for given parameters (beta, J, B, c).
    print("The values of the message 'u' and magnetization 'm_0' are determined by solving the following self-consistency equations:")
    print("m_0 = tanh(beta * (B + c * u))")
    print("tanh(beta * u) = tanh(beta * J) * tanh(beta * (B + (c - 1) * u))")

print_susceptibility_formula()