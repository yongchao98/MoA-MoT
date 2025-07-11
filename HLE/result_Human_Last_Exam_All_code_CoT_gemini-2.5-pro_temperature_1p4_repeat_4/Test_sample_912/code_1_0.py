def solve_work_cycle():
    """
    This function prints the derived formula for the work done by the current source.
    The formula is derived from the principles of electromagnetism and thermodynamics
    for the given electromechanical system and cycle.
    """
    
    # The derived formula for the work done W by the current source per cycle is:
    # W = - ( (mu - mu_0) * N^2 * w * (x2 - x1) * (I2^2 - I1^2) ) / (2 * g)
    
    # We print the formula part by part as requested.
    # Note: mu represents μ, mu_0 represents μ₀, N, w, g, x1, x2, I1, I2 are the given parameters.
    
    term1 = "-(mu - mu_0)"
    term2 = "N^2"
    term3 = "w"
    term4 = "(x2 - x1)"
    term5 = "(I2^2 - I1^2)"
    denominator = "2*g"

    print("The work done by the current source for each cycle is:")
    print(f"W = ({term1} * {term2} * {term3} * {term4} * {term5}) / ({denominator})")
    
solve_work_cycle()