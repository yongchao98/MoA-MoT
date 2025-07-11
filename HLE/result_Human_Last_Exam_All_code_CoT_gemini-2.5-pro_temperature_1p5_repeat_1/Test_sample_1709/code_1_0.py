import numpy as np

def solve(n, k):
    """
    Calculates the number of non-zero elements of the vector field V(P).
    
    Based on the analysis, the complex definitions simplify significantly under a set of reasonable assumptions,
    leading to a simple formula for the number of non-zero elements.
    
    The logic is as follows:
    1. The integrals defining the matrices are assumed to evaluate to 1, -1, and 0 to ensure P lies on the manifold.
    2. The integral equation for the vector field V simplifies, revealing V as the gradient of a function F.
    3. For the case n=k, the number of non-zero elements of V(P) is n^2 / 2.
    4. This result is generalized to the formula nk/2 for n >= k.
    
    Args:
        n (int): An even integer >= 6.
        k (int): An even integer >= 6, with n >= k.
    
    Returns:
        int: The number of non-zero elements, l(n,k).
    """
    
    if not (isinstance(n, int) and isinstance(k, int)):
        raise TypeError("n and k must be integers.")
    if not (n >= 5 and k >= 5):
        raise ValueError("n and k must be greater than or equal to 5.")
    if not (n % 2 == 0 and k % 2 == 0):
        raise ValueError("n and k must be even numbers.")
    if not n >= k:
        raise ValueError("n must be greater than or equal to k.")

    # The number of non-zero elements is given by the formula nk/2.
    l_nk = n * k // 2
    
    print(f"Given n = {n} and k = {k}")
    print(f"The number of non-zero elements l(n,k) is calculated by the formula n * k / 2.")
    print(f"l({n},{k}) = {n} * {k} / 2 = {l_nk}")
    
# Example usage with some even numbers n, k >= 6 and n >= k.
# You can change these values to test other cases.
n_example = 10
k_example = 6

solve(n_example, k_example)
