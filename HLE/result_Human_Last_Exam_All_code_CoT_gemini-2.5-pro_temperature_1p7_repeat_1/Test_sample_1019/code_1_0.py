import math

def calculate_hamiltonicity_threshold(n, eta):
    """
    Calculates the d-threshold for Hamiltonicity for a graph H_n U G(n,p).

    Args:
        n (int): The number of vertices in the graph.
        eta (float): The deficit parameter, where d = n/2 - eta.
    
    Returns:
        float: The threshold probability p.
    """
    # Check if the inputs are valid according to the problem description
    if not isinstance(n, int) or n <= 0:
        print("Error: n must be a positive integer.")
        return None
    if not (0.5 <= eta <= n / 64):
        print(f"Error: eta must be in the range [0.5, n/64]. For n={n}, this range is [0.5, {n/64}].")
        print(f"The provided value eta={eta} is outside this range.")
        return None

    # The formula for the threshold p
    # p = 2 * eta / (n * (n - 2 * eta - 1))
    
    numerator = 2 * eta
    denominator_part = n - 2 * eta - 1
    
    if denominator_part <= 0:
        print("Error: The denominator becomes non-positive with the given n and eta.")
        return None
        
    denominator = n * denominator_part
    p = numerator / denominator

    # Output the result showing each number in the equation
    print("Based on the known sharp threshold for this problem, the formula for p is:")
    print("p = (2 * eta) / (n * (n - 2 * eta - 1))\n")
    print(f"For the given values n = {n} and eta = {eta}:")
    
    print(f"p = (2 * {eta}) / ({n} * ({n} - 2 * {eta} - 1))")
    print(f"p = {numerator} / ({n} * {denominator_part})")
    print(f"p = {numerator} / {denominator}")
    print(f"p = {p:.10f}")
    
    return p

if __name__ == '__main__':
    # Example values for n and eta satisfying the constraints.
    # Let's choose a large n and an eta within its valid range.
    n_example = 2000
    eta_example = 20  # This satisfies 0.5 <= 20 <= 2000/64 = 31.25

    print(f"--- Calculating Hamiltonicity Threshold for n={n_example}, eta={eta_example} ---")
    calculate_hamiltonicity_threshold(n_example, eta_example)
