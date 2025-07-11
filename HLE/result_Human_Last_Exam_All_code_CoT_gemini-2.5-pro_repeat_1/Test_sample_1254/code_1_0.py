import math

def calculate_H(a, b, c, d, r, t):
    """
    Calculates the upper bound H based on the derived formula.
    
    H(a, b, c, d, r, t) = (|a| * b * t) / (c * d**2 * r)
    
    Args:
        a (float): The parameter k (k < 0).
        b (float): The L1 norm of rho, ||rho(0,.)||_L1.
        c (float): The value of pi.
        d (float): The cutoff radius nu.
        r (float): The value of the function rho(tau, x).
        t (float): The time upper limit of integration.
        
    Returns:
        float: The calculated value of the upper bound H.
    """
    if a > 0:
        print("Warning: Parameter 'a' (k) is expected to be negative.")
    if d <= 0:
        raise ValueError("Parameter 'd' (nu) must be positive.")
    if r <= 0:
        raise ValueError("Parameter 'r' (rho) must be positive.")
    
    abs_a = abs(a)
    d_squared = d**2
    numerator = abs_a * b * t
    denominator = c * d_squared * r
    
    H = numerator / denominator
    
    print("Calculating the upper bound H:")
    print(f"H = (|a| * b * t) / (c * d^2 * r)")
    print(f"H = (|{a}| * {b} * {t}) / ({c} * {d}^2 * {r})")
    print(f"H = ({abs_a} * {b} * {t}) / ({c} * {d_squared} * {r})")
    print(f"H = {numerator} / {denominator}")
    print(f"H = {H}")
    
    return H

if __name__ == '__main__':
    # Example values for the parameters
    # a = k
    k = -2.0
    # b = ||rho(0,.)||_L1
    rho_L1_norm = 1.0
    # c = pi
    pi_val = math.pi
    # d = nu
    nu = 0.1
    # r = rho(tau, x)
    rho_val = 0.5
    # t = t
    time = 10.0

    print("Given parameters:")
    print(f"a = k = {k}")
    print(f"b = ||rho(0,.)||_L1 = {rho_L1_norm}")
    print(f"c = pi = {pi_val}")
    print(f"d = nu = {nu}")
    print(f"r = rho(tau, x) = {rho_val}")
    print(f"t = {time}\n")

    # Calculate and print the result
    H_value = calculate_H(k, rho_L1_norm, pi_val, nu, rho_val, time)
