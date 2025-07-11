import math

def get_cost(X):
    """
    Finds the lexicographically smallest non-negative integer pair (a, b)
    such that (1-a)(1-b) = X and a, b != 1.
    This function uses a directed search for efficiency.
    """
    best_pair = None
    # We search for pairs (a,b) that result in X.
    # A small search radius is sufficient because we are looking for the lexicographically smallest pair.
    search_limit = abs(X) + 5
    for a in range(search_limit):
        if a == 1:
            continue
        # We need (1-a)*(1-b) = X, so 1-b = X/(1-a)
        if (1 - a) != 0 and X % (1 - a) == 0:
            one_minus_b = X // (1 - a)
            b = 1 - one_minus_b
            if b >= 0 and b != 1:
                current_pair = (a, b)
                if best_tuple is None or current_pair < best_pair:
                    best_pair = current_pair
    return best_pair


# Based on the analysis, the minimal number of manifolds is l=3.
# The condition is sum(X_i) = (l-1)/2 = 1.
l = 3
target_sum = 1
best_tuple = None

# We search for a partition of the target_sum into l non-zero integers {X_1, ..., X_l}
# that minimizes the lexicographical cost.
# A bounded search is sufficient to find the optimal solution.
search_range = 5
for x1 in range(-search_range, search_range + 1):
    if x1 == 0: continue
    for x2 in range(-search_range, search_range + 1):
        if x2 == 0: continue
        
        x3 = target_sum - x1 - x2
        if x3 == 0: continue
        
        partition = [x1, x2, x3]
        
        costs = []
        possible = True
        for x in partition:
            cost = get_cost(x)
            if cost is None:
                possible = False
                break
            costs.append(cost)
        
        if not possible:
            continue
            
        # Sort the pairs to form the lexicographically smallest final tuple
        costs.sort()
        
        candidate_tuple = tuple(item for pair in costs for item in pair)
        
        if best_tuple is None or candidate_tuple < best_tuple:
            best_tuple = candidate_tuple

# The discovered best tuple is used to formulate the final answer and verification.
a1, b1 = best_tuple[0], best_tuple[1]
a2, b2 = best_tuple[2], best_tuple[3]
a3, b3 = best_tuple[4], best_tuple[5]

chi1 = 4 * (1 - a1) * (1 - b1)
chi2 = 4 * (1 - a2) * (1 - b2)
chi3 = 4 * (1 - a3) * (1 - b3)
chi_final = chi1 + chi2 + chi3 - 2 * (l - 1)

# Print the verification equation as requested.
print(f"The equation for the Euler characteristic of the connect-sum is:")
print(f"{chi1} + {chi2} + ({chi3}) - 2*({l}-1) = {chi_final}")

# Print the final tuple in the requested format.
print("The lexicographically least tuple is:")
print(f"({','.join(map(str, best_tuple))})")