import math

def calculate_simplices(N, k, max_n):
    """
    Calculates and prints the number of n-simplices for n from 0 to max_n.

    The number of n-simplices in N_dot(Z_N)_k/ is the number of non-decreasing
    sequences of n+1 integers in the range [k, N]. This is a combination with
    repetition problem.

    Args:
        N (int): The upper bound of the integer range.
        k (int): The starting object for the over category.
        max_n (int): The maximum dimension of the simplices to calculate.
    """
    # The number of available integers to choose from is m.
    m = N - k + 1
    
    print(f"Given N = {N}, k = {k}.")
    print(f"The number of objects to choose from is m = N - k + 1 = {N} - {k} + 1 = {m}.")
    print("The number of n-simplices is given by the formula for combinations with repetition:")
    print("C(m + r - 1, r), where r = n + 1.")
    print(f"This simplifies to C({m} + n, n + 1).\n")

    results = []
    for n in range(max_n + 1):
        # r is the number of items to choose in the combination formula, which is n+1.
        r = n + 1
        
        # C_n_val and C_k_val are the parameters for the binomial coefficient C(n, k).
        C_n_val = m + r - 1
        C_k_val = r
        
        # Calculate the result using math.comb for accuracy.
        num_simplices = math.comb(C_n_val, C_k_val)
        results.append(num_simplices)

        # Build the strings for the expanded formula representation.
        numerator_terms = [str(C_n_val - i) for i in range(C_k_val)]
        denominator_terms = [str(C_k_val - i) for i in range(C_k_val)]

        numerator_str = " * ".join(numerator_terms)
        denominator_str = " * ".join(denominator_terms)
        
        print(f"For n = {n}:")
        if len(numerator_terms) > 1:
            print(f"  Number of simplices = C({C_n_val}, {C_k_val}) = ({numerator_str}) / ({denominator_str}) = {num_simplices}")
        else:
            print(f"  Number of simplices = C({C_n_val}, {C_k_val}) = {numerator_str} / {denominator_str} = {num_simplices}")

    return results

# Problem parameters
N_val = 200
k_val = 13
n_val = 5

final_answers = calculate_simplices(N_val, k_val, n_val)
# The final answer is returned after all print statements.
# This format is requested for the final result.
# print(f"\n<<<{final_answers}>>>")