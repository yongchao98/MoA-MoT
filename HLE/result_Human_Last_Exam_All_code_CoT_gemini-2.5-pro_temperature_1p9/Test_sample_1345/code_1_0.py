import math

def calculate_max_zeros(N):
    """
    Calculates the maximal number of complex zeros for the given matrix B(k).
    
    Args:
        N: The dimension of the matrix.
        
    Returns:
        The maximal number of complex zeros.
    """
    if not isinstance(N, int) or N < 1:
        print("Error: N must be a positive integer.")
        return None
        
    if N == 1:
        result = 0
        print(f"For N = {N}, the polynomial for k1 has degree 1, so its root must be real.")
        print(f"The number of complex zeros is 0.")
    else:
        # For N >= 2, the maximal number is given by the formula N * 2**(N-1)
        base = 2
        exponent = N - 1
        # It's good practice to use math.pow for potentially large exponents,
        # but ** is fine for typical integer values.
        term2 = base**exponent
        result = N * term2
        
        print(f"For N = {N}, the maximal number of complex zeros is given by the formula N * 2^(N-1).")
        print("The calculation is:")
        print(f"{N} * {base}^({N} - 1) = {N} * {base}^{exponent} = {N} * {term2} = {result}")

    return result

# --- Main execution ---
# Let's calculate the result for a specific N, for instance N=4, as an example.
N_example = 4
final_answer = calculate_max_zeros(N_example)
# The final answer in the requested format is the numerical result for N=4.
# If you wish to test other integer values for N, you can change N_example.

if final_answer is not None:
    # Per instructions, returning the final answer in the specified format.
    # The printed output above explains the general formula and the derivation.
    # This final line provides the specific answer for the chosen N_example.
    pass # The thought process will conclude with the specific format string.
