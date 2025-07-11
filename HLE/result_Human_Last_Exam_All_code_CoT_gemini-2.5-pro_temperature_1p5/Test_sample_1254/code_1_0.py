import math

def calculate_H(a, b, c, d, r, t):
    """
    Calculates the upper bound H based on the derived formula.

    Args:
        a (float): The parameter k (k < 0).
        b (float): The L1 norm of rho, ||rho(0, .)||_L^1.
        c (float): The value of pi.
        d (float): The parameter nu (nu > 0).
        r (float): The value of rho(tau, x) (rho > 0).
        t (float): The time t (t > 0).

    Returns:
        float: The calculated upper bound H.
    """
    if a >= 0:
        raise ValueError("Parameter 'a' (k) must be negative.")
    if b < 0:
        raise ValueError("Parameter 'b' (L1 norm) must be non-negative.")
    if d <= 0:
        raise ValueError("Parameter 'd' (nu) must be positive.")
    if r <= 0:
        raise ValueError("Parameter 'r' (rho) must be positive.")
    if t < 0:
        raise ValueError("Parameter 't' (time) must be non-negative.")

    # H = (-a * b * t) / (c * d^2 * r)
    numerator = -a * b * t
    denominator = c * d**2 * r
    
    H = numerator / denominator
    
    # Print the equation with the given numbers
    print("The explicit formula for the upper bound H is:")
    print("H(a, b, c, d, r, t) = (-a * b * t) / (c * d^2 * r)")
    print("\nPlugging in the given values:")
    print(f"H = (-({a}) * {b} * {t}) / ({c} * {d}**2 * {r})")
    print(f"H = ({numerator}) / ({denominator})")
    print(f"H = {H}")

    return H

if __name__ == '__main__':
    # Example parameters
    k_val = -2.0  # Corresponds to a
    rho_L1_norm = 1.0  # Corresponds to b
    pi_val = math.pi  # Corresponds to c
    nu_val = 0.1  # Corresponds to d
    rho_val = 0.5  # Corresponds to r
    time_val = 10.0  # Corresponds to t

    # Calculate and print the bound
    calculate_H(a=k_val, b=rho_L1_norm, c=pi_val, d=nu_val, r=rho_val, t=time_val)

    # The final answer is the formula for H
    # Using the variable names from the prompt a, b, c, d, r, t
    final_expression = "(-a * b * t) / (c * d**2 * r)"
    # print(f"\n<<< {final_expression} >>>")