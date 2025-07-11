import numpy as np

def calculate_log_density_term(d, lam, x_dot_p):
    """
    Calculates a single term of the log-likelihood ratio, log(h(v(x))).
    
    Args:
        d (int): The dimension.
        lam (float): The lambda parameter.
        x_dot_p (float): The dot product between x and p.
        
    Returns:
        float: The value of log(h(v(x))).
    """
    # t = ||v|| = arccos(x.p)
    # Clamp the argument to arccos to avoid domain errors due to floating point inaccuracies
    t = np.arccos(np.clip(x_dot_p, -1.0, 1.0))
    
    # sinc(t) = sin(t)/t. For t=0, the limit is 1.
    if t == 0:
        log_sinc_t = 0
    else:
        log_sinc_t = np.log(np.sin(t) / t)
        
    log_h = -t**2 / (2 * lam) + (d - 2) * log_sinc_t
    return log_h, t

def solve_log_ratio(d, lam):
    """
    Calculates and prints the value of l(d, lambda).
    """
    if not isinstance(d, int) or d < 4:
        raise ValueError("d must be an integer greater than or equal to 4.")
    if not isinstance(lam, (int, float)) or lam < 1:
        raise ValueError("lambda must be a number greater than or equal to 1.")

    # Dot products for x1 and x2 with p = 1_d/sqrt(d)
    x1_dot_p = np.sqrt(3 / d)
    x2_dot_p = np.sqrt(2 / d)

    # Calculate log density terms and the norms t1, t2
    log_h1, t1 = calculate_log_density_term(d, lam, x1_dot_p)
    log_h2, t2 = calculate_log_density_term(d, lam, x2_dot_p)

    # Calculate the final log ratio
    ell = log_h1 - log_h2

    # Output the numbers in the final equation as requested
    print(f"Calculating l(d, lambda) for d={d}, lambda={lam}")
    print("-" * 30)
    print(f"Equation: l(d, lambda) = log(h(x1)) - log(h(x2))")
    print(f"where log(h(x)) = -t^2/(2*lambda) + (d-2)*log(sinc(t)) and t = arccos(p.x)")
    print("-" * 30)
    print(f"For x1: p.x1 = sqrt(3/{d}) = {x1_dot_p:.4f}")
    print(f"t1 = arccos({x1_dot_p:.4f}) = {t1:.4f}")
    print(f"log(h(x1)) = -({t1:.4f})^2/(2*{lam}) + ({d}-2)*log(sinc({t1:.4f})) = {log_h1:.4f}")
    print("-" * 30)
    print(f"For x2: p.x2 = sqrt(2/{d}) = {x2_dot_p:.4f}")
    print(f"t2 = arccos({x2_dot_p:.4f}) = {t2:.4f}")
    print(f"log(h(x2)) = -({t2:.4f})^2/(2*{lam}) + ({d}-2)*log(sinc({t2:.4f})) = {log_h2:.4f}")
    print("-" * 30)
    print(f"Final Result: l({d}, {lam}) = {log_h1:.4f} - ({log_h2:.4f}) = {ell}")
    
    return ell

if __name__ == '__main__':
    # Using d=4 and lambda=1 as example values from the function's domain.
    d_val = 4
    lambda_val = 1.0
    final_answer = solve_log_ratio(d_val, lambda_val)
    # The final answer is also printed inside the function.
    # To conform to the output format, we capture and print it again here if needed,
    # or just let the function handle all printing.
