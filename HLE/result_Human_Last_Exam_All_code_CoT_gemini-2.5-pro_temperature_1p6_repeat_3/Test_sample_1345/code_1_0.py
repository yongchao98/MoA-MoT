import math

def calculate_max_complex_zeros(N):
    """
    Calculates the maximal possible number of complex zeros for the given matrix B(k).

    The number of zeros is given by a polynomial in k_1 of degree d.
    The degree d is N * 2^(N-1).
    For N=1, the degree is 1, so the single root must be real. Thus, 0 complex zeros.
    For N>=2, the degree is even, and for a generic matrix A, all roots can be complex.
    So the maximal number of complex zeros is the degree of the polynomial.
    
    Args:
        N (int): The dimension of the matrix. Must be a positive integer.
        
    Returns:
        int: The maximal number of complex zeros.
    """
    if not isinstance(N, int) or N <= 0:
        raise ValueError("N must be a positive integer.")

    print(f"For N = {N}:")
    
    if N == 1:
        result = 0
        print("The degree of the polynomial is 1, which implies the only root is real.")
        print(f"The maximal number of complex zeros is {result}.")
    else:
        # The formula for the maximal number of complex zeros is N * 2^(N-1) for N>=2
        term1 = N
        term2 = 2
        term3 = N - 1
        
        # In the final code you still need to output each number in the final equation!
        print("The calculation for the maximal number of complex zeros follows the formula: N * 2^(N-1)")
        print("Substituting the values:")
        
        result = term1 * (term2 ** term3)
        
        print(f"{term1} * ({term2} ** {term3}) = {result}")
        
    return result

# Example: Calculate the number for a specific N.
# You can change the value of N here.
try:
    N = 3
    calculate_max_complex_zeros(N)
except ValueError as e:
    print(e)
