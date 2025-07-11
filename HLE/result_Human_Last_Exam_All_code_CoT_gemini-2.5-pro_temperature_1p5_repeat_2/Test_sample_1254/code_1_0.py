import math

def solve_H(a, b, c, d, r, t):
    """
    Calculates the upper bound H.
    
    Args:
        a: a float representing k (<0)
        b: a float representing the L1 norm of rho
        c: a float representing pi
        d: a float representing nu (>0)
        r: a float representing the infimum of rho(tau, x) (>0)
        t: a float representing the time t (>0)
    """

    if a >= 0:
        raise ValueError("Parameter 'a' (k) must be negative.")
    if b < 0:
        raise ValueError("Parameter 'b' (L1 norm of rho) must be non-negative.")
    if d <= 0:
        raise ValueError("Parameter 'd' (nu) must be positive.")
    if r <= 0:
        raise ValueError("Parameter 'r' (infimum of rho) must be positive.")
    if t < 0:
        raise ValueError("Parameter 't' must be non-negative.")
        
    numerator = -a * b * t
    denominator = c * d**2 * r
    
    H = numerator / denominator
    
    # Output the final equation with numerical values
    print("The explicit expression for the upper bound H is: H = (-a * b * t) / (c * d**2 * r)")
    print(f"H = (-({a}) * {b} * {t}) / ({c} * {d}**2 * {r})")
    print(f"H = {numerator} / {denominator}")
    print(f"H = {H}")

# Example usage with placeholder values.
# The user can replace these values with their specific ones.
k_val = -2.0
rho_L1_norm = 1.0
pi_val = math.pi
nu_val = 0.1
rho_infimum = 0.5
t_val = 10.0

solve_H(a=k_val, b=rho_L1_norm, c=pi_val, d=nu_val, r=rho_infimum, t=t_val)

# The final answer format as requested by the user prompt.
final_expression = "(-a * b * t) / (c * d**2 * r)"
print(f"\n<<<H(a, b, c, d, r, t) = {final_expression}>>>")
