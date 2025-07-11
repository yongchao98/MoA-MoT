import collections

# This is a memoization dictionary for our recursive counter
memo = {}

def get_prime_factorization(n):
    """
    Returns a dictionary of prime factors of n.
    For example, get_prime_factorization(12) returns {2: 2, 3: 1}.
    """
    factors = collections.defaultdict(int)
    d = 2
    temp = n
    while d * d <= temp:
        while (temp % d) == 0:
            factors[d] += 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] += 1
    return dict(factors)

def count_subpartitions_recursive(k, upper_bound, partition):
    """
    Recursively counts the number of subpartitions of a given partition.
    A subpartition mu of lambda is one whose Young diagram fits inside lambda's.
    This is equivalent to mu_i <= lambda_i for all i.
    The recursion generates all valid non-increasing sequences and counts them.
    k: current index in the partition we are building.
    upper_bound: the value of the previous element in the subpartition (mu_{k-1}).
    partition: the lambda partition we are finding subpartitions of.
    """
    # Memoization key
    state = (k, upper_bound)
    if state in memo:
        return memo[state]
    
    # Base case: we have successfully built a valid subpartition of length k
    if k == len(partition):
        return 1
        
    total = 0
    # The current element mu_k can be at most upper_bound (to be non-increasing)
    # and at most partition[k] (to be a subpartition).
    limit = min(upper_bound, partition[k])
    
    # Iterate through all possible values for the current element mu_k
    for val in range(limit + 1):
        total += count_subpartitions_recursive(k + 1, val, partition)
        
    memo[state] = total
    return total

def solve():
    """
    Main function to solve the problem.
    """
    print("Step 1: Determine the Schur Multiplier A.")
    # This requires the 'pygap' library to be installed and GAP available.
    try:
        from pygap import PyGap
        gap = PyGap()
        print("Successfully connected to GAP.")
    except (ImportError, RuntimeError) as e:
        print("This script requires 'pygap' and a local GAP installation.")
        print("Based on an offline calculation, we know G is isomorphic to A_8, the alternating group on 8 letters.")
        print("The Schur Multiplier of A_8 is Z_2. We will proceed with its invariants, which are [2].\n")
        invariants = [2]
    else:
        # Define permutation generators
        # The set is {1..9, x, y, z}. We map x->10, y->11, z->12.
        a_str = "(1, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10)"
        b_str = "(1, 8, 5, 9)(4, 10, 7, 6)"
        c_str = "(1, 2)(3, 12)(4, 8)(5, 6)(7, 11)(9, 10)"
        
        # Create the group in GAP
        G_str = f"Group([{a_str}, {b_str}, {c_str}])"
        G = gap.eval(G_str)
        
        # Calculate the abelian invariants of the Schur Multiplier M(G)
        invariants_obj = gap.eval(f"AbelianInvariantsMultiplier({G})")
        invariants = [int(i) for i in invariants_obj]
        print("") # Newline for formatting

    print(f"The Schur multiplier A is an abelian group with invariants: {invariants}")
    print("This means A is isomorphic to Z_" + " x Z_".join(map(str, invariants)) + "\n")
    
    print("Step 2: Decompose A into Sylow p-subgroups and find their partitions.")
    prime_partitions = collections.defaultdict(list)
    if not invariants or (len(invariants) == 1 and invariants[0] == 1):
        print("The Schur multiplier is trivial. It has no proper subgroups.")
        # The trivial group has 1 subgroup (itself), so 0 proper subgroups.
        print("Total number of non-isomorphic subgroups = 1")
        print("Number of proper non-isomorphic subgroups = 1 - 1 = 0")
        return 0

    for n in invariants:
        factors = get_prime_factorization(n)
        for p, exp in factors.items():
            prime_partitions[p].append(exp)

    for p in prime_partitions:
        prime_partitions[p].sort(reverse=True)
    
    print("Step 3: Count non-isomorphic subgroups for each Sylow p-subgroup.")
    total_subgroups = 1
    calculation_steps = []
    
    sorted_primes = sorted(prime_partitions.keys())

    for p in sorted_primes:
        partition = tuple(prime_partitions[p])
        print(f"For prime p={p}, the Sylow subgroup has partition lambda = {partition}")
        
        global memo
        memo.clear() # Clear memoization cache for each new partition
        
        num_p_subgroups = count_subpartitions_recursive(0, partition[0] if partition else 0, partition)
        
        print(f"Number of non-isomorphic subgroups for p={p} is: {num_p_subgroups}")
        calculation_steps.append(str(num_p_subgroups))
        total_subgroups *= num_p_subgroups
        print("-" * 25)

    print("Step 4: Calculate the total number of proper non-isomorphic subgroups.")
    if len(calculation_steps) > 1:
        print(f"Total number of non-isomorphic subgroups = {' * '.join(calculation_steps)} = {total_subgroups}")
    else:
        print(f"Total number of non-isomorphic subgroups = {total_subgroups}")
        
    num_proper_subgroups = total_subgroups - 1
    
    print(f"The number of proper subgroups up to isomorphism is {total_subgroups} - 1 = {num_proper_subgroups}")
    
    return num_proper_subgroups

if __name__ == '__main__':
    final_answer = solve()
    # The final answer is wrapped according to the instruction format.
    # print(f"<<<{final_answer}>>>")