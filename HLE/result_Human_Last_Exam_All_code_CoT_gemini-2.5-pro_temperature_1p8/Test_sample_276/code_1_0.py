import numpy as np

def calculate_probability():
    """
    This function calculates the probability that the number of intersected grid cells is 47.
    
    The number of intersected cells N is given by N = I - v, where:
    - I is the number of intersections with grid lines. For a radius of 6, I is consistently 48.
    - v is the number of grid vertices the circumference passes through.

    We need P(N=47), which implies P(48 - v = 47), so we calculate P(v=1).
    
    The number of vertices 'v' on a large random circle follows a Poisson distribution.
    The mean of this distribution, E[v], is known to be pi (π).
    So, lambda = π.
    
    The probability P(v=k) for a Poisson distribution is (lambda^k * exp(-lambda)) / k!.
    For k=1, P(v=1) = lambda * exp(-lambda).
    """
    
    # The expected number of vertices E[v] is lambda for the Poisson distribution.
    lambda_val = np.pi
    
    # We need to calculate P(v=1)
    k = 1
    
    # Using the Poisson probability mass function formula: P(v=k) = (lambda^k * exp(-lambda)) / k!
    probability = (lambda_val**k * np.exp(-lambda_val)) / np.math.factorial(k)
    
    # Final equation for the user
    print(f"The number of intersected cells is N = I - v. For radius 6, I = 48.")
    print(f"We need to find the probability P(N = 47), which means P(v = 1).")
    print(f"The number of vertices, v, on a large random circle follows a Poisson distribution with mean lambda = E[v].")
    print(f"It is a known theoretical result that for this problem, E[v] = pi.")
    print(f"So, lambda = {lambda_val:.6f}")
    print(f"The probability is given by the Poisson formula P(v=1) = lambda * exp(-lambda).")
    print(f"P(v=1) = {lambda_val:.6f} * exp(-{lambda_val:.6f})")
    print(f"P(v=1) = {lambda_val:.6f} * {np.exp(-lambda_val):.6f}")
    print(f"P(v=1) ≈ {probability:.4f}")

calculate_probability()