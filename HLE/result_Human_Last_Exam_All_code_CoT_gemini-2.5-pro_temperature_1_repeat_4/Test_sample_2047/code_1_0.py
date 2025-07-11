import numpy as np

def solve_l_d_lambda(d, lam):
    """
    This function calculates the value of l(d, lambda) based on the derived formula.

    The formula is: l(d, lambda) = ln [ (g(pi/2) + g(a1)) / (g(pi/2) + g(a2)) ]
    where:
    g(a) = exp(-a^2 / (2*lambda)) * (sinc(a))^(d-2)
    a1 = arccos(sqrt(3/d))
    a2 = arccos(sqrt(2/d))
    """

    # Helper function for sinc(x) = sin(x)/x
    def sinc(x):
        if np.isclose(x, 0.0):
            return 1.0
        else:
            return np.sin(x) / x

    # Helper function for g(a), the proportional density based on norm 'a'
    def g(a, d_val, lam_val):
        # Clamp sinc_val to avoid log(0) or negative bases for the power
        sinc_val = sinc(a)
        if sinc_val <= 0:
            # This case happens for a > pi, but our norms are in [0, pi/2].
            # Added for robustness.
            return 0.0
        
        term1 = np.exp(-a**2 / (2 * lam_val))
        term2 = sinc_val**(d_val - 2)
        return term1 * term2

    print(f"Calculating l(d, lambda) for d = {d} and lambda = {lam}\n")

    # Step 1: Calculate the pre-image norms a1 and a2
    a1 = np.arccos(np.sqrt(3 / d))
    a2 = np.arccos(np.sqrt(2 / d))
    a_pi_half = np.pi / 2

    # Step 2: Calculate the density contributions g(a) for each norm
    g_pi_half = g(a_pi_half, d, lam)
    g_a1 = g(a1, d, lam)
    g_a2 = g(a2, d, lam)

    # Step 3: Sum contributions to get the proportional densities for x1 and x2
    p_x1_prop = g_pi_half + g_a1
    p_x2_prop = g_pi_half + g_a2
    
    # Check for non-positive values before taking log
    if p_x1_prop <= 0 or p_x2_prop <= 0:
        print("Error: Proportional density is not positive. Cannot compute log-ratio.")
        return None
        
    # Step 4: Calculate the ratio of densities
    ratio = p_x1_prop / p_x2_prop

    # Step 5: Calculate the final value l(d, lambda) = ln(ratio)
    l_value = np.log(ratio)
    
    # Output the numbers in the final equation as requested
    print("--- Final Equation Components ---")
    print(f"a1 = arccos(sqrt(3/{d})) = {a1:.6f}")
    print(f"a2 = arccos(sqrt(2/{d})) = {a2:.6f}")
    print(f"a_pi_half = pi/2 = {a_pi_half:.6f}")
    print(f"g(a1) = {g_a1:.6f}")
    print(f"g(a2) = {g_a2:.6f}")
    print(f"g(pi/2) = {g_pi_half:.6f}")
    print(f"p_prop(x1) = g(pi/2) + g(a1) = {p_x1_prop:.6f}")
    print(f"p_prop(x2) = g(pi/2) + g(a2) = {p_x2_prop:.6f}")
    print(f"Ratio = p_prop(x1) / p_prop(x2) = {ratio:.6f}")
    print("---------------------------------")
    print(f"Final result: l({d}, {lam}) = ln(Ratio) = {l_value:.6f}")

    return l_value

if __name__ == '__main__':
    # As d and lambda are not specified in the prompt, we use the smallest valid integer values:
    # d must be >= 4
    # lambda must be >= 1
    d_val = 4
    lambda_val = 1
    
    final_answer = solve_l_d_lambda(d_val, lambda_val)
    # The final answer in the required format
    if final_answer is not None:
        print(f"<<<{final_answer:.7f}>>>")
