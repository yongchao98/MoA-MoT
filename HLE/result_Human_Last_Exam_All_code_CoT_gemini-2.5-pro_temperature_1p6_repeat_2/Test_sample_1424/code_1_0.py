def print_green_function_equation():
    """
    This function prints the equation for the bare Green's function G_0
    as a function of frequency (ω) and single-particle energy (ϵ_k).
    The equation is constructed and printed part by part.
    """
    
    # Using Unicode for better readability of Greek letters
    omega = "ω"
    epsilon_k = "ϵₖ"
    mu = "µ"
    delta = "δ"
    theta = "θ"

    # Print the equation part by part to meet the output requirements
    # "G₀(ω, ϵₖ) = "
    print("G₀({0}, {1}) = ".format(omega, epsilon_k), end="")

    # First term: Particle propagator (for states above the Fermi level)
    # "θ(ϵₖ - µ)"
    print("{0}({1} - {2})".format(theta, epsilon_k, mu), end="")
    # " / "
    print(" / ", end="")
    # "(ω - ϵₖ + iδ)"
    print("({0} - {1} + i{2})".format(omega, epsilon_k, delta), end="")

    # Separator
    # " + "
    print(" + ", end="")

    # Second term: Hole propagator (for states below the Fermi level)
    # "θ(µ - ϵₖ)"
    print("{0}({1} - {2})".format(theta, mu, epsilon_k), end="")
    # " / "
    print(" / ", end="")
    # "(ω - ϵₖ - iδ)"
    print("({0} - {1} - i{2})".format(omega, epsilon_k, delta))

print_green_function_equation()