def solve_set_problem():
    """
    This function demonstrates the construction of an infinite set 'x' that is
    almost disjoint from every set in a given countable collection 'S'.
    """

    # Step 1: Define a sample countable collection S of infinite sets.
    # For this example, let S be composed of sets of multiples.
    # s_k = { k, 2k, 3k, 4k, ... } for k = 1, 2, 3, ...
    # We will work with a finite portion of this infinite collection.
    num_sets_in_S = 50
    S = []
    # To ensure sets are 'infinite enough' for the demo, we make them large.
    max_val_in_sets = 100000
    for i in range(1, num_sets_in_S + 1):
        # We use a list to keep elements sorted, which is like the
        # enumerating function f_k(n).
        s_i = sorted(list({j * i for j in range(1, max_val_in_sets // i)}))
        if s_i: # Ensure the set is not empty
            S.append(s_i)

    # Let's verify the setup
    print(f"Demonstrating the construction with a collection S of {len(S)} infinite sets.")
    print(f"s_0 = {S[0][:10]}...")
    print(f"s_1 = {S[1][:10]}...")
    print("-" * 20)

    # Step 2: Implement the construction of the set x = {a_0, a_1, ...}
    # We will construct the first `n_elements` of x.
    n_elements = 25
    x = []
    a_prev = -1  # Represents a_{n-1}

    print(f"Constructing the first {n_elements} elements of the set x...")

    # The loop from n=0 to n_elements-1 calculates a_n
    for n in range(n_elements):
        # The number of sets we diagonalize against at step n is n+1 (s_0 to s_n)
        num_sets_to_consider = n + 1
        if num_sets_to_consider > len(S):
            print(f"Stopping at n={n} as we've run out of sample sets in S.")
            break
            
        # Get the n-th element from each set s_0, ..., s_n. This is f_k(n).
        # We need to handle the case where n is out of bounds for a list.
        # This can happen if our example sets aren't large enough.
        elements_at_n = []
        valid_construction = True
        for k in range(num_sets_to_consider):
            if n >= len(S[k]):
                print(f"Warning: set s_{k} does not have an element at index {n}. This demo requires larger sets.")
                valid_construction = False
                break
            elements_at_n.append(S[k][n])
        
        if not valid_construction:
            break
            
        # bound(n) = max(f_0(n), ..., f_n(n))
        bound_n = max(elements_at_n)

        # a_n = max(a_{n-1}, bound(n)) + 1
        a_n = max(a_prev, bound_n) + 1
        x.append(a_n)
        a_prev = a_n

        # The following line implements the full definition with each number printed out
        # This might get a bit verbose for larger n.
        if n < 5: # Only show the details for the first few steps
          print(f"Step n={n}:")
          print(f"  a_{n-1} = {x[n-1] if n>0 else -1}")
          print(f"  Elements at index n={n} from s_0 to s_{n}: {elements_at_n}")
          print(f"  bound({n}) = max({elements_at_n}) = {bound_n}")
          print(f"  a_{n} = max(a_{n-1}, bound({n})) + 1 = max({x[n-1] if n>0 else -1}, {bound_n}) + 1 = {a_n}")
          print("")


    print("-" * 20)
    print("Construction complete.")
    print(f"The first {len(x)} elements of the constructed set x are:")
    print(x)

solve_set_problem()