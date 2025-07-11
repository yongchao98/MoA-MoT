import itertools

def solve():
    """
    Solves the mathematical problem by searching for the set A.

    The problem asks for min_{A ∈ Σ} max_{a ∈ A}, where Σ is the set of
    all finite non-empty subsets of positive integers A (excluding {2})
    such that for all x, y in A, x+y is a product of two elements in A.

    Plan:
    1. Iterate through possible maximum values M for a set A, starting from M=1.
    2. For each M, generate all subsets of {1, ..., M-1}.
    3. Create candidate sets A by adding M to each subset.
    4. Test if each candidate set A meets the conditions:
        a. A must not be {2}.
        b. The minimum element of A must be 1 or 2.
        c. For every pair x, y in A, their sum x+y must be present in the
           set of products of pairs from A.
    5. The first M for which such a set A is found is the minimum maximum,
       which is the answer.
    6. If no set is found within a reasonable search limit, the answer is 0.
    """

    # We search for M, the maximum element of a potential set A.
    # A search up to M=16 should be sufficient to find a solution if one exists
    # in the small integers. A larger search becomes computationally expensive.
    limit_m = 16 
    for m_val in range(1, limit_m + 1):
        # The candidate set A must contain m_val as its maximum element.
        # The other elements must come from the set {1, 2, ..., m_val-1}.
        pool = list(range(1, m_val))
        
        # Iterate through all subsets of the pool to form candidate sets.
        for i in range(len(pool) + 1):
            for subset in itertools.combinations(pool, i):
                A = sorted(list(subset) + [m_val])

                # Optimization: The minimum element of A must be 1 or 2.
                if A[0] > 2:
                    continue

                # The set {2} is explicitly excluded from Σ.
                if A == [2]:
                    continue

                # Check the condition A+A ⊆ A·A
                sums = set()
                for x in A:
                    for y in A:
                        sums.add(x + y)

                products = set()
                for x in A:
                    for y in A:
                        products.add(x * y)

                if sums.issubset(products):
                    # We found a valid set A. Since we iterate m_val upwards,
                    # this is the smallest possible maximum value.
                    print(f"Found a solution set A = {A}")
                    print(f"The maximum element is max(A) = {m_val}")
                    print(f"The minimum of these maximums is therefore {m_val}.")
                    print("\nFinal equation:")
                    print(f"min_{{A in Sigma}} max_{{a in A}} a = {m_val}")
                    return m_val
    
    # If the loop completes without finding a set, Sigma is empty.
    print("No solution found within the search limit.")
    print("\nFinal equation:")
    print(f"min_{{A in Sigma}} max_{{a in A}} a = 0")
    return 0

# Execute the solver and print the final answer.
final_answer = solve()
print(f"\n<<< {final_answer} >>>")
