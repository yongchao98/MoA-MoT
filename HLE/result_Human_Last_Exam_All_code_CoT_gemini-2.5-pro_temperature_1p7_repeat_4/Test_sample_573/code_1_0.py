import math

def calculate_simplices(N, k, n_max):
    """
    Calculates and prints the number of n-simplices for n from 0 to n_max.

    Args:
        N (int): The maximum value for the subcategory Z_N.
        k (int): The object over which the slice category is taken.
        n_max (int): The maximum dimension of the simplices to calculate.
    """
    print(f"Calculating the number of n-simplices for N={N}, k={k}, and n up to {n_max}.")
    print("The number of n-simplices is given by the formula: C(N - k + n + 1, n + 1)")
    print("-" * 60)
    
    results = []
    for n in range(n_max + 1):
        # The general formula is C(N - k + n + 1, n + 1)
        comb_n = N - k + n + 1
        comb_k = n + 1
        
        # Calculate the result using math.comb
        try:
            result = math.comb(comb_n, comb_k)
            results.append(result)
            
            # Print the equation and the result
            print(f"For n = {n}:")
            print(f"Number of simplices = C({N} - {k} + {n} + 1, {n} + 1) = C({comb_n}, {comb_k}) = {result}")
        except ValueError:
            print(f"Could not calculate for n={n}, C({comb_n}, {comb_k}) is not valid.")

    print("-" * 60)
    print("Final list of results for n=0 to 5:")
    print(results)


# Given parameters
N_val = 200
k_val = 13
n_max_val = 5

calculate_simplices(N_val, k_val, n_max_val)
