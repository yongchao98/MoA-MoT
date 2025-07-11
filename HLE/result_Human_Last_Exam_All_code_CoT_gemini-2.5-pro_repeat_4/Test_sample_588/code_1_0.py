import collections

def find_minimal_solution():
    """
    Finds the lexicographically least tuple (a_1, b_1, ..., a_l, b_l)
    with l minimal, such that M(a_i, b_i) is not full, but their
    connect-sum is full.
    """
    # l must be odd and greater than 1. Minimal l is 3.
    ell = 3
    target_x_sum = (ell - 1) / 2

    # Generate a list of possible pairs (a,b) and their x values,
    # sorted lexicographically. a, b are non-negative, not 1.
    # We enforce a <= b for a canonical representation of pairs.
    limit = 5 # A search limit sufficient for this problem
    candidate_pairs = []
    for a in range(limit):
        if a == 1:
            continue
        for b in range(a, limit):
            if b == 1:
                continue
            x = (1 - a) * (1 - b)
            candidate_pairs.append(((a, b), x))

    # Systematically search for the lexicographically smallest solution.
    # The first one found by this ordered search will be the minimal one.
    for i in range(len(candidate_pairs)):
        p1, x1 = candidate_pairs[i]
        for j in range(i, len(candidate_pairs)):  # p2 >= p1
            p2, x2 = candidate_pairs[j]
            
            required_x3 = target_x_sum - x1 - x2
            
            # Find smallest p3 >= p2 that satisfies the equation
            for k in range(j, len(candidate_pairs)): # p3 >= p2
                p3, x3 = candidate_pairs[k]
                if x3 == required_x3:
                    # Solution found.
                    solution_pairs = sorted([p1, p2, p3])
                    
                    # Print the verification of the equation
                    x_vals = [(1-p[0])*(1-p[1]) for p in solution_pairs]
                    print(f"Minimal l = {ell}")
                    print("The equation is: 2 * sum((1-a_i)*(1-b_i)) = l - 1")
                    print(f"2 * ( (1-{solution_pairs[0][0]})*(1-{solution_pairs[0][1]}) + (1-{solution_pairs[1][0]})*(1-{solution_pairs[1][1]}) + (1-{solution_pairs[2][0]})*(1-{solution_pairs[2][1]}) ) = {ell} - 1")
                    print(f"2 * ( {x_vals[0]} + {x_vals[1]} + {x_vals[2]} ) = {ell - 1}")
                    print(f"2 * ( {sum(x_vals)} ) = {ell - 1}")
                    print(f"{2 * sum(x_vals)} = {ell - 1}")
                    print("-" * 20)

                    # Unpack into the final flat tuple
                    result_tuple = tuple(item for pair in solution_pairs for item in pair)
                    
                    # Print the final answer as a flat tuple
                    final_answer_str = "(" + ",".join(map(str, result_tuple)) + ")"
                    print("Lexicographically least tuple:")
                    print(final_answer_str)
                    return

find_minimal_solution()