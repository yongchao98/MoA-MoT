import itertools

def solve_set_problem():
    """
    Solves the user's problem by searching for the required set.

    The problem asks to evaluate min(max(A)) over a set Σ. The definition
    of Σ involves the condition A+A ⊆ A × A. As this is not possible
    under a literal interpretation (integers vs. ordered pairs), we
    interpret A × A as the set of outputs from a pairing function applied
    to pairs of elements from A. We use the Cantor pairing function for this.

    The code will search for the smallest 'm' such that there exists a set
    A ∈ Σ with max(A) = m.
    """
    
    # Using Cantor pairing function for positive integers x, y ≥ 1
    # π(x, y) = 1/2 * (x + y - 2)(x + y - 1) + y
    def pairing_function(x, y):
        return (x + y - 2) * (x + y - 1) // 2 + y

    def check_set(A):
        """
        Checks if a given set A is in Σ.
        A ∈ Σ if:
        1. A is a finite subset of positive integers (implicit in representation).
        2. A is not ∅ or {2}.
        3. A+A ⊆ π(A,A).
        """
        # Condition 2
        if not A or A == {2}:
            return False, None, None

        # Condition 3
        sum_set = {a + b for a in A for b in A}
        pair_set = {pairing_function(a, b) for a in A for b in A}

        if sum_set.issubset(pair_set):
            return True, sum_set, pair_set
        else:
            return False, sum_set, pair_set

    # Analysis shows min(A) must be 1. So we search for sets A of the
    # form {1} U S U {m} where S is a subset of {2, ..., m-1}.
    m = 1
    while True:
        middle_elements = list(range(2, m))
        
        # Iterate through all subsets of {2, ..., m-1}
        for i in range(len(middle_elements) + 1):
            for subset in itertools.combinations(middle_elements, i):
                # Construct candidate set A
                A = {1, m} | set(subset)
                
                # Check if it works
                is_in_sigma, sum_set, pair_set = check_set(A)
                
                if is_in_sigma:
                    print(f"The set Σ is not empty.")
                    print(f"Found a qualifying set A = {sorted(list(A))}.")
                    print(f"The maximum element in this set is {m}.")
                    print(f"Since we searched in increasing order of the maximum element, this is the minimum possible maximum.")
                    print("\n--- Condition Verification ---")
                    # Final equation details requested by user prompt
                    print(f"A+A = {sorted(list(sum_set))}")
                    print(f"π(A,A) = {sorted(list(pair_set))}")
                    print(f"Check: A+A ⊆ π(A,A) is True.")
                    
                    result = m
                    print(f"\nThe result min_{{A ∈ Σ}} max_{{a ∈ A}} a is {result}.")
                    return

        m += 1
        if m > 20: # Safety break
            print("Search terminated without finding a solution.")
            return

# Run the solver
solve_set_problem()