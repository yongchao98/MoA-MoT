def print_overlap_integral_formula():
    """
    Prints the analytical expression for the overlap integral S(2s, 2s)
    for two hydrogenic 2s orbitals.
    """
    
    zeta_R = "\u03B6R" # Unicode for zeta followed by R, representing rho
    
    print("The analytical expression for the overlap integral S between two 2s orbitals is a function of the internuclear distance R and the effective nuclear charge \u03B6.")
    print(f"Defining a dimensionless variable \u03C1 = {zeta_R}, the formula is:")
    print()
    
    # The derived formula S(rho) = exp(-rho/2) * (1 + rho/2 + rho^2/12 + rho^4/240)
    # The print statement below constructs this formula string, showing each coefficient.
    
    formula = f"S(\u03C1) = exp(-\u03C1/2) * (1 + (1/2)\u03C1 + (1/12)\u03C1\u00B2 + (1/240)\u03C1\u2074)"
    
    print(formula)
    print()
    print("where:")
    print("  S is the overlap integral")
    print("  \u03B6 (zeta) is the effective nuclear charge")
    print("  R is the internuclear distance")
    print("  exp() is the exponential function")
    print("  \u00B2 and \u2074 represent the superscripts for squared and to the power of 4")

print_overlap_integral_formula()