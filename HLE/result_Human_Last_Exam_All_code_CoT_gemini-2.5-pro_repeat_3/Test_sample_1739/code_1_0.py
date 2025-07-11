import math

def solve_bubble_frequency():
    """
    Calculates the 3rd term of the nonlinear frequency correction for the Rayleigh-Plesset equation.
    """
    # Step 1: Assume a value for the polytropic index gamma.
    # A common value for air bubbles used in such problems is 4/3.
    gamma = 4/3

    # Step 2: Calculate the linear oscillation frequency w0.
    # w0^2 = 3 * gamma
    w0 = math.sqrt(3 * gamma)

    # Step 3: The nonlinear frequency correction coefficient (w2/a^2) has three terms.
    # Term 1 = w0^5 / 6
    # Term 2 = 3 * w0^3 / 6
    # Term 3 = 3 * w0 / 6
    
    term1_val = (w0**5) / 6
    term2_val = (3 * w0**3) / 6
    term3_val = (3 * w0) / 6
    
    total_val = term1_val + term2_val + term3_val

    # Step 4: Output the final equation with evaluated numbers.
    print("The equation for the nonlinear frequency correction coefficient is:")
    print(f"{term1_val:.4f} + {term2_val:.4f} + {term3_val:.4f} = {total_val:.4f}")
    
    # Step 5: Output the value of the 3rd term.
    print("\nThe 3rd term of the nonlinear correction to the linear oscillation frequency is:")
    print(int(term3_val))

solve_bubble_frequency()