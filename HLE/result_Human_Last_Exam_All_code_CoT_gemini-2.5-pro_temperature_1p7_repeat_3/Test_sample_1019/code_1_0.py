import math

def calculate_hamiltonicity_threshold(n, eta):
    """
    Calculates the d-threshold for Hamiltonicity based on the worst-case
    graph construction of two disjoint cliques.

    Args:
        n (int): The number of vertices.
        eta (float): The parameter defining the minimum degree d = n/2 - eta.
    
    Returns:
        float: The threshold probability p.
    """
    
    # Check if eta is in the given range for validity
    if not (0.5 <= eta <= n / 64):
        print("Warning: eta is outside the specified range [1/2, n/64].")

    # Sizes of the two disjoint cliques in the worst-case construction
    size_A = n / 2 - eta + 1
    size_B = n / 2 + eta - 1

    # The threshold is determined by the need to connect the smaller clique (A)
    # to the larger one (B). The critical probability is approximately ln(|A|) / |B|.
    if size_A <= 1:
        # log is not defined for <=0, and threshold would be different for size_A=1
        raise ValueError("n and eta values result in a trivial clique size.")

    p_threshold = math.log(size_A) / size_B
    
    # Output the equation with the specific numbers
    print(f"The calculation for n={n} and eta={eta} is:")
    print(f"p = ln({size_A}) / {size_B}")
    print(f"p = {math.log(size_A):.5f} / {size_B}")
    print(f"p = {p_threshold}")
    
    return p_threshold

# Example values
n_val = 20000
eta_val = 100

# Calculate and print the result. The final value will be captured for the answer.
result = calculate_hamiltonicity_threshold(n_val, eta_val)

# The final answer is the numeric result of the calculation.
# print(f'<<<{result:.5f}>>>')