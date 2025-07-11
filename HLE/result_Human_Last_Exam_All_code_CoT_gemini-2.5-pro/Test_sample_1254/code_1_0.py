import math

def calculate_H_bound(a, b, c, d, r, t):
    """
    Calculates the upper bound H based on the derived formula.

    The problem asks to estimate the upper bound H for the quantity:
    | integral from 0 to t of (f(tau, x) / rho(tau, x)) d(tau) |

    Under the simplifying assumption that rho is time-independent, i.e., rho(tau, x) = rho(x),
    the upper bound H is given by the formula:
    H(a, b, c, d, r, t) = (-a * b * t) / (c * d**2 * r)

    Args:
        a (float): The parameter k, which must be negative (a < 0).
        b (float): The L1 norm of rho, ||rho(0, .)||_L1, which must be non-negative.
        c (float): The value of pi.
        d (float): The parameter nu, the radius of the excluded ball, which must be positive (d > 0).
        r (float): The value of rho(x), which must be positive (r > 0).
        t (float): The time duration, which must be non-negative.
    
    Returns:
        float: The calculated value of the upper bound H.
    """
    if a >= 0:
        raise ValueError("Parameter 'a' (k) must be negative.")
    if d <= 0:
        raise ValueError("Parameter 'd' (nu) must be positive.")
    if r <= 0:
        raise ValueError("Parameter 'r' (rho(x)) must be positive.")
    if b < 0:
        raise ValueError("Parameter 'b' (L1 norm of rho) must be non-negative.")
    if t < 0:
        raise ValueError("Parameter 't' (time) must be non-negative.")
        
    numerator = -a * b * t
    denominator = c * (d**2) * r
    
    if denominator == 0:
        raise ValueError("Denominator cannot be zero.")
        
    H = numerator / denominator
    
    print("This script calculates the upper bound H based on the provided parameters.")
    print("The formula for H is: H = (-a * b * t) / (c * d^2 * r)")
    print("-" * 30)
    print(f"Given parameters:")
    print(f"  a = k = {a}")
    print(f"  b = ||rho(0, .)||_L1 = {b}")
    print(f"  c = pi = {c}")
    print(f"  d = nu = {d}")
    print(f"  r = rho(x) = {r}")
    print(f"  t = {t}")
    print("-" * 30)
    
    print("Calculation steps:")
    print(f"  Numerator = -a * b * t = -({a}) * {b} * {t} = {numerator}")
    print(f"  Denominator = c * d^2 * r = {c} * {d}^2 * {r} = {denominator}")
    print(f"  H = Numerator / Denominator = {numerator} / {denominator}")
    print("-" * 30)
    print(f"Final Result: H = {H}")
    
    return H

# --- Example Usage ---
# You can change these values to see how the bound changes.
k_val = -2.0
rho_L1_norm_val = 1.0
pi_val = math.pi
nu_val = 0.1
rho_x_val = 0.5
t_val = 10.0

calculate_H_bound(k_val, rho_L1_norm_val, pi_val, nu_val, rho_x_val, t_val)