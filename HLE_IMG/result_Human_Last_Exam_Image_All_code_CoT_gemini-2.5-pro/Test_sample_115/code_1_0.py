from itertools import product

def get_possible_p_values():
    """Generates all possible values for P(A), which is a product of 4 integers from {1, 2, 3}."""
    possible_values = set()
    for combo in product([1, 2, 3], repeat=4):
        p = combo[0] * combo[1] * combo[2] * combo[3]
        possible_values.add(p)
    return sorted(list(possible_values))

def find_min_max_models(n_is):
    """
    Finds the min and max number of models (k) such that n_is can be written
    as a sum of k valid P values.
    """
    p_values = get_possible_p_values()
    
    # --- Find Maximum Models ---
    # To maximize k, we use the smallest possible P value, which is 1.
    max_k = n_is 

    # --- Find Minimum Models ---
    # Find the smallest k for which a solution exists.
    min_k = float('inf')
    
    # We can use dynamic programming to find the minimum number of terms.
    # dp[i] will store the minimum number of P values that sum to i.
    dp = [float('inf')] * (n_is + 1)
    dp[0] = 0
    
    for i in range(1, n_is + 1):
        for p in p_values:
            if i >= p:
                if dp[i - p] != float('inf'):
                    dp[i] = min(dp[i], dp[i - p] + 1)
    
    min_k = dp[n_is]

    return (min_k, max_k)

# The number of independent sets of size 4 in the given graph is 7.
N_IS = 7

min_models, max_models = find_min_max_models(N_IS)

print(f"The number of independent sets of size 4 is N_IS = {N_IS}.")
print(f"We need to find min/max k such that N_IS = P_1 + P_2 + ... + P_k.")
print(f"Maximum number of models: {max_models} (by using P_i = 1 for all i)")
print(f"Minimum number of models: {min_models} (e.g., 7 can be written as 6 + 1, where 6=3*2*1*1 and 1=1*1*1*1 are valid P values)")
print(f"The pair (min, max) is ({min_models}, {max_models}).")
