import collections

def get_prime_factor_sets(n_max):
    """
    Computes the set of prime factors for each number from 1 to n_max.
    """
    factor_sets = collections.defaultdict(set)
    for i in range(2, n_max + 1):
        if not factor_sets[i]:  # i is prime
            for j in range(i, n_max + 1, i):
                factor_sets[j].add(i)
    # Add empty set for n=1 explicitly if needed, but defaultdict(set) handles it.
    return factor_sets

def count_allowed_pairs(n_max):
    """
    Counts the number of allowed pairs (a,b) with 1 <= a, b <= n_max.
    """
    # Step 1: Get prime factor sets for all numbers up to n_max
    prime_factors = get_prime_factor_sets(n_max)
    prime_factors[1] = set()

    # Step 2: Group numbers by their prime factor sets and count them
    groups = collections.defaultdict(int)
    for n in range(1, n_max + 1):
        # Use frozenset as dict key
        p_set = frozenset(prime_factors[n])
        groups[p_set] += 1

    # Step 3: Get the list of unique prime factor sets and their counts
    unique_sets = list(groups.keys())
    m = len(unique_sets)
    
    # Step 4: Calculate the total number of allowed pairs
    allowed_pairs_count = 0
    for i in range(m):
        for j in range(m):
            set_i = unique_sets[i]
            set_j = unique_sets[j]
            
            # Check the condition for an allowed pair
            if set_i.issubset(set_j) or set_j.issubset(set_i):
                count_i = groups[set_i]
                count_j = groups[set_j]
                allowed_pairs_count += count_i * count_j
                
    return allowed_pairs_count

# Set the limit for a and b
N = 1000

# Calculate the result
result = count_allowed_pairs(N)
print(result)
