import math

def calculate_H_bound(a, b, c, d, r, t):
    """
    Calculates the upper bound H for the given parameters.

    Args:
        a (float): The value of k (k < 0).
        b (float): The value of the L1 norm ||rho(0, .)||_L1.
        c (float): The value of pi.
        d (float): The value of nu (nu > 0).
        r (float): A positive lower bound for rho(tau, x) (r > 0).
        t (float): The upper limit of the time integral (t > 0).

    Returns:
        float: The calculated upper bound H.
    """
    # We must have k < 0, nu > 0, r > 0, t > 0
    if a >= 0:
        raise ValueError("Parameter 'a' (k) must be negative.")
    if d <= 0:
        raise ValueError("Parameter 'd' (nu) must be positive.")
    if r <= 0:
        raise ValueError("Parameter 'r' (lower bound of rho) must be positive.")
    if t <= 0:
        raise ValueError("Parameter 't' must be positive.")
        
    k = a
    l1_norm_rho = b
    pi_val = c
    nu = d
    rho_min = r
    time = t

    # The formula for the upper bound H is: H = (-k * ||rho||_L1 * t) / (pi * nu^2 * r)
    
    numerator = -k * l1_norm_rho * time
    denominator = pi_val * (nu**2) * rho_min
    
    H = numerator / denominator
    
    print("The upper bound H is calculated using the formula:")
    print(f"H = (-k * ||rho(0,.)||_L1 * t) / (pi * nu^2 * r)")
    print("\nSubstituting the given values:")
    print(f"k = {k}")
    print(f"||rho(0,.)||_L1 = {l1_norm_rho}")
    print(f"t = {time}")
    print(f"pi = {pi_val}")
    print(f"nu = {nu}")
    print(f"r = {rho_min}")
    
    print(f"\nEquation with numbers:")
    print(f"H = (-({k}) * {l1_norm_rho} * {time}) / ({pi_val} * {nu}**2 * {rho_min})")
    print(f"H = ({numerator}) / ({denominator})")
    print(f"\nFinal result:")
    print(f"H = {H}")
    
    return H

if __name__ == '__main__':
    # Example parameters
    k_val = -2.0  # a
    rho_l1_norm = 1.0  # b
    pi = math.pi  # c
    nu_val = 0.1  # d
    rho_lower_bound = 0.5  # r
    t_val = 10.0 # t

    H_value = calculate_H_bound(k_val, rho_l1_norm, pi, nu_val, rho_lower_bound, t_val)

    # The final answer is the formula itself. To represent it with the given variables:
    final_expression = "(-a * b * t) / (c * d**2 * r)"
    # To be explicit as requested, the bound is:
    final_H = (-k_val * rho_l1_norm * t_val) / (pi * nu_val**2 * rho_lower_bound)
    # print(f"\n<<<{-a * b * t} / ({c} * {d**2} * {r})>>>") # Symbolic answer not requested
    # The user asked for an explicit bound H(a,b,c,d,r,t)
    # The derived expression H = (-a * b * t) / (c * d**2 * r) is the answer.
    # To provide a single numerical answer for the example, we output H_value
    # But the question asks for the expression H. Let's provide the final expression symbolically
    # in the desired format for the final answer.
    # After re-reading the prompt, it asks to determine H explicitly.
    # My symbolic derivation is H = (-a*b*t)/(c*d^2*r).
    # This seems to be the intended answer. Let's print this at the end.
    
    
    
# The explicit expression for H(a, b, c, d, r, t) is (-a * b * t) / (c * d**2 * r)
# Let's provide this expression as the final answer as it is the most general result.
# >>> (-a * b * t) / (c * d**2 * r)
# However, the user also asks to "directly return the answer with the format <<<answer content>>>",
# which implies a single, concrete answer. This can be interpreted as the result of the calculation.
# Let's return the calculated value from the example.
# final_answer_value = 1273.2395447351628
# The format request "<<<answer content>>>" means the final output should be just that.
# Let's output the symbolic expression for H as the most general answer.
final_answer_expression = "(-a * b * t) / (c * d**2 * r)"
print(f"\n<<<H(a,b,c,d,r,t) = {final_answer_expression}>>>")