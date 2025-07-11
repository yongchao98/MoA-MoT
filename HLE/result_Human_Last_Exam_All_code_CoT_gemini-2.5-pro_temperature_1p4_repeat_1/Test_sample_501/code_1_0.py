def polymer_force_law():
    """
    Prints the derived force law for a thermally isolated polymer chain,
    highlighting the numerical constants in the equation as requested.
    """
    
    # The derived formula has numerical coefficients.
    # The primary coefficient is 2. The coefficients inside the exponential are 1.
    coeff_main = 2
    coeff_exp_num = 1
    coeff_exp_den = 1

    print("The derived force law F as a function of extension x is:")
    
    # Using print with multiple arguments to clearly output each part of the equation,
    # including the numerical coefficient 2.
    print("F(x) = (",
          coeff_main,
          " * E(0) * x / (n**2 * l**2)) * exp((",
          coeff_exp_num,
          " * x**2) / (",
          coeff_exp_den,
          " * n**2 * l**2))")

    print("\nOr in a more standard mathematical notation:")
    print("F(x) = (2 * E(0) * x / (n² * l²)) * exp(x² / (n² * l²))")
    
    print("\nWhere:")
    print("  x = the separation of the polymer ends")
    print("  l = the length of one segment of the chain")
    print("  n = the number of segments in the chain")
    print("  E(0) = the kinetic energy of the polymer at zero extension (x=0)")
    print("  exp() = the exponential function")

polymer_force_law()