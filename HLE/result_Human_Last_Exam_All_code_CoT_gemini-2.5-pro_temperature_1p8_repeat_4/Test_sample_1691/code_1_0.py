import math

def solve():
    """
    This function calculates the coefficients for the analytical approximation of the integral I(epsilon)
    and prints the final formula.
    """
    # The integral is I(epsilon) = integral from 0 to 15 of 1 / (epsilon + 9*x^5 + 5*x^6 + 9*x^8) dx.
    # For small epsilon, the behavior is dominated by the term with the lowest power of x in the denominator,
    # which is 9.0 * x^5.0.
    
    # Let's denote the leading term as a * x^n.
    a = 9.0
    n = 5.0
    
    # The asymptotic formula for an integral of the form integral[1 / (epsilon + a*x^n)] dx for small epsilon
    # is given by I(epsilon) approx C * epsilon^p.
    
    # The power 'p' is derived from the scaling property and is equal to -(n-1)/n.
    p = -(n - 1.0) / n
    
    # The coefficient 'C' is given by C = a^(-1/n) * integral from 0 to infinity of 1/(1+u^n) du.
    # The standard integral is integral[1/(1+u^n)]du from 0 to inf = (pi/n) / sin(pi/n).
    # So, C = a^(-1/n) * (pi/n) / sin(pi/n).
    
    # Calculate the coefficient C
    try:
        coefficient_c = (a**(-1.0/n)) * (math.pi / n) / math.sin(math.pi / n)
    except (ValueError, ZeroDivisionError) as e:
        print(f"Error calculating coefficient C: {e}")
        return

    # Print the derived formula and the values of C and p.
    print("The analytical formula that approximates I(epsilon) for the small epsilon regime is of the form:")
    print("I(epsilon) \u2248 C * \u03B5^p")
    print("\nBased on the leading term 9.0 * x^5.0, the constants are:")
    
    print(f"\nThe coefficient C is calculated as:")
    print(f"C = {a}^(-1/{n}) * (pi/{n}) / sin(pi/{n})")
    print(f"C \u2248 {coefficient_c}")

    print(f"\nThe power p is calculated as:")
    print(f"p = -({n}-1)/{n}")
    print(f"p = {p}")
    
    # Output the final equation with numerical values.
    print("\nThus, the final analytical approximation is:")
    print(f"I(\u03B5) \u2248 {coefficient_c} * \u03B5^({p})")
    
    # Outputting the coefficient C for the final answer block.
    print(f"\n<<<{coefficient_c}>>>")

solve()