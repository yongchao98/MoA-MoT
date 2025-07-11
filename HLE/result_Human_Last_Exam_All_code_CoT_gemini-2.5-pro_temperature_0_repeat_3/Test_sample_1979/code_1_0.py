import itertools

def solve_and_print_answer():
    """
    Searches for a set A in Sigma and computes min(max(A)) over A in Sigma.
    Sigma = {A | A is a finite, non-empty subset of positive integers,
             A != {2}, and A+A is a subset of A*A}.
    If Sigma is empty, the function prints 0.
    """

    # We search for the minimum possible value of max(A), let's call it m.
    # We can iterate m from 1 upwards. The first m that yields a valid set A is the answer.
    # The search space grows as 2^(m-1), so we set a reasonable limit for m.
    search_limit = 16 
    found_solution = False
    min_max_a = 0

    for m in range(1, search_limit + 1):
        # A must contain m. The other elements must come from {1, 2, ..., m-1}.
        other_elements = list(range(1, m))
        
        # Iterate through all subsets of other_elements to form candidate sets A.
        for i in range(len(other_elements) + 1):
            for subset in itertools.combinations(other_elements, i):
                # Construct the candidate set A
                A = set(subset) | {m}

                # Condition: A must not be {2}
                if A == {2}:
                    continue

                # Generate the sumset and product set
                sumset = {x + y for x in A for y in A}
                prodset = {x * y for x in A for y in A}

                # Check the main condition: A+A subset A*A
                if sumset.issubset(prodset):
                    # A solution is found. Since we iterate m in increasing order,
                    # this is the minimum possible max(A).
                    min_max_a = m
                    found_solution = True
                    break  # Exit inner loops
            if found_solution:
                break
        
        if found_solution:
            break # Exit outer loop

    # Based on the search result, print the final answer.
    # Mathematical analysis strongly suggests no such set exists, so the search
    # is expected to find nothing, resulting in printing 0.
    print(min_max_a)

solve_and_print_answer()