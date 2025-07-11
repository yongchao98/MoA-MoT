import numpy as np

def dimensionless_lindhard_function(x):
    """
    Calculates the dimensionless part of the static Lindhard function.
    The argument x is the ratio q / (2*k_F).
    """
    # Avoid division by zero at x=0 by handling the limit case.
    if np.abs(x) < 1e-12:
        return 1.0
    # Avoid issues with the logarithm for x >= 1
    if np.abs(x) >= 1:
        # For x > 1, the argument of the log is negative. We use the absolute value as per the formula.
        # For x = 1, there's a singularity, but we are interested in the q->0 limit.
        x_val = np.abs(x)
    else:
        x_val = x
    
    # Using np.log1p(y) is more numerically stable than np.log(1+y) for small y
    log_term = (np.log1p(x_val) - np.log1p(-x_val))
    
    return 0.5 + ((1 - x_val**2) / (4 * x_val)) * log_term

def solve_lindhard_problem():
    """
    Solves for the numerical value of the Lindhard function at q=0.
    """
    # The Lindhard function at (q=0, w=0) is -g(eps_F) * F(0), where F is the
    # dimensionless part of the function. We calculate F(0) by taking the limit.
    # We use a very small value for x to approximate the limit x -> 0.
    x_limit = 1e-9
    limit_value_F = dimensionless_lindhard_function(x_limit)

    # The result for Pi(0,0) is -1 * g(eps_F) * F(0).
    # The universal numerical factor is -1 * F(0).
    numerical_factor = -1.0
    
    final_value = numerical_factor * limit_value_F

    print("The static Lindhard function Pi(q, 0) can be written as:")
    print("Pi(q, 0) = -g(eps_F) * F(q / 2*k_F)")
    print("where g(eps_F) is the density of states at the Fermi level.")
    print("\nWe need to evaluate this at q=0, which means we need the limit of F(x) as x->0.")
    print(f"\nCalculating the limit of F(x) as x -> 0 (using x = {x_limit}):")
    print(f"F(0) = {limit_value_F:.6f}")

    print("\nThe Lindhard function at q=0 is therefore:")
    print("Pi(0, 0) = -1 * g(eps_F) * F(0)")
    
    print("\nThe numerical value of the ratio Pi(0, 0) / g(eps_F) is -1 * F(0).")
    print("\nFinal Equation:")
    print(f"{numerical_factor} * {limit_value_F:.1f} = {final_value:.1f}")
    
    print(f"\nThus, the Lindhard polarization function at zero frequency and zero momentum transfer, expressed in units of the density of states at the Fermi level g(eps_F), is -1.")


solve_lindhard_problem()