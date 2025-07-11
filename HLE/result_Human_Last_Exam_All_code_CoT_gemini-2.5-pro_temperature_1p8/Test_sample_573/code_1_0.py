import math

def calculate_simplices(N, k):
    """
    Calculates the number of n-simplices for the over-category N.(Z_N)_k/
    for n from 0 to 5.
    
    An n-simplex in this category corresponds to a sequence (j_0, ..., j_n)
    such that k <= j_0 <= j_1 <= ... <= j_n <= N.
    
    The number of such sequences is given by the combination with repetition formula:
    C((N - k + 1) + (n + 1) - 1, n + 1) = C(N - k + n + 1, n + 1)
    """
    
    print(f"Calculating the number of n-simplices for N = {N} and k = {k}:\n")
    
    results = []
    for n in range(6):
        # Upper argument of the binomial coefficient
        upper_arg = N - k + n + 1
        # Lower argument of the binomial coefficient
        lower_arg = n + 1
        
        # Calculate the number of simplices
        num_simplices = math.comb(upper_arg, lower_arg)
        results.append(num_simplices)
        
        print(f"For n = {n}:")
        print(f"Number of {n}-simplices = C({N} - {k} + {n} + 1, {n} + 1) = C({upper_arg}, {lower_arg}) = {num_simplices}")
        print("-" * 20)
        
    return results

if __name__ == "__main__":
    N = 200
    k = 13
    
    # Run the calculation and store the results
    simplex_counts = calculate_simplices(N, k)
    
    # The final answer is requested in a specific format
    # In this case, we'll represent it as a list of counts for n=0 through n=5.
    final_answer = simplex_counts
    # The submission format below is for the final, single answer.
    # print(f"\nFinal Answer: {final_answer}")