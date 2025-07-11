import math

def calculate_lower_bound(N, q):
    """
    Calculates the lower bound for m based on the derived formula.

    Args:
        N (int): The number of rows in the input matrix X.
        q (int): The sparsity parameter.
    """
    if not isinstance(N, int) or not isinstance(q, int) or N <= 0 or q <= 0:
        print("Error: N and q must be positive integers.")
        return
    if q > N / 2:
        print(f"Warning: The derivation assumes q <= N/2. The provided q={q} does not satisfy this.")
    
    # The lower bound is floor(N / (2q))
    m_lower_bound = math.floor(N / (2.0 * q))
    
    print(f"Given N = {N} and q = {q}:")
    print(f"The number of independent 'tasks' we can construct is M = floor(N / (2q)).")
    print(f"M = floor({N} / (2 * {q})) = floor({N / (2*q):.2f}) = {int(m_lower_bound)}")
    print(f"The dimension of the hidden layer 'm' must be at least M.")
    print(f"Thus, the lower bound for m is {int(m_lower_bound)}.")
    print(f"Asymptotically, the lower bound for m is Omega(N/q).")

# Example usage with some hypothetical values for N and q
# Note: The user did not provide specific values, so we use examples.
print("--- Example 1 ---")
calculate_lower_bound(N=1000, q=10)
print("\n--- Example 2 ---")
calculate_lower_bound(N=500, q=5)
print("\n--- Example 3 ---")
calculate_lower_bound(N=200, q=1)
