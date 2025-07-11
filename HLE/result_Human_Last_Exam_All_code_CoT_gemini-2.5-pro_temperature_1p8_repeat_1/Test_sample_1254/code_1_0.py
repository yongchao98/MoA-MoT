import math

def calculate_H_bound(k, rho_L1_norm, pi_val, nu, rho_min, t):
    """
    Calculates the explicit upper bound H based on the derived formula.

    The formula for the upper bound H is:
    H = (-k * ||rho(0,.)||_L1 * t) / (pi * nu^2 * r)

    Args:
        k (float): The parameter k, which must be negative.
        rho_L1_norm (float): The L1 norm of the initial function rho, ||rho(0,.)||_L1.
        pi_val (float): The value of pi.
        nu (float): The radius nu for the modified Riesz transform, must be positive.
        rho_min (float): A positive lower bound for the function rho(tau, x).
        t (float): The time t, must be non-negative.

    Returns:
        float: The calculated value of the upper bound H.
    """
    # Validate the input parameters based on the problem's constraints.
    if k >= 0:
        raise ValueError("Parameter 'k' must be negative.")
    if rho_L1_norm < 0:
        raise ValueError("Parameter 'rho_L1_norm' (b) must be non-negative.")
    if nu <= 0:
        raise ValueError("Parameter 'nu' (d) must be positive.")
    if rho_min <= 0:
        raise ValueError("Parameter 'rho_min' (r) must be positive.")
    if t < 0:
        raise ValueError("Parameter 't' must be non-negative.")

    numerator = -k * rho_L1_norm * t
    denominator = pi_val * (nu ** 2) * rho_min

    return numerator / denominator

def main():
    """
    Main function to demonstrate the calculation of H.
    """
    # Set example values for the parameters, corresponding to H(a,b,c,d,r,t)
    # a = k_val
    k_val = -2.5
    # b = rho_L1_norm_val
    rho_L1_norm_val = 1.0
    # c = pi_val
    pi_val = math.pi
    # d = nu_val
    nu_val = 0.1
    # r = rho_min_val
    rho_min_val = 0.05
    # t = t_val
    t_val = 10.0

    print("--- Calculating the upper bound H ---")
    print("The formula for H is: H = (-k * ||rho(0,.)||_L1 * t) / (pi * nu^2 * r)")
    
    # Outputting each number in the final equation
    print("\nUsing the following parameter values:")
    print(f"k = {k_val}")
    print(f"||rho(0,.)||_L1 = {rho_L1_norm_val}")
    print(f"t = {t_val}")
    print(f"pi = {pi_val}")
    print(f"nu = {nu_val}")
    print(f"r = {rho_min_val}")

    # Display the equation with the numerical values substituted
    print("\nThe final equation with numbers is:")
    # We use .4f for pi to keep it clean, but use full precision in calculation.
    print(f"H = (-({k_val}) * {rho_L1_norm_val} * {t_val}) / ({pi_val:.4f} * {nu_val}^2 * {rho_min_val})")
    
    try:
        # Calculate the result
        H_value = calculate_H_bound(k_val, rho_L1_norm_val, pi_val, nu_val, rho_min_val, t_val)
        print("\n--- Result ---")
        print(f"The calculated upper bound H is: {H_value}")
    except ValueError as e:
        print(f"\nError in calculation: {e}")

if __name__ == "__main__":
    main()
