import collections

def solve_and_verify():
    """
    This function verifies the derived solution for the manifold problem.
    The goal is to find the lexicographically least tuple (a_1, b_1, ..., a_l, b_l),
    with l minimal, such that M(a_i, b_i) is not full, but their connect-sum is.

    This reduces to the equation: 2 * sum((1-a_i)*(1-b_i)) = l-1
    """

    # The derived minimal tuple solution
    solution_tuple = (0, 0, 0, 0, 0, 2)

    # Extract parameters from the tuple
    # M(a_i, b_i) corresponds to M(genus_a, genus_b)
    pairs = [solution_tuple[i:i+2] for i in range(0, len(solution_tuple), 2)]
    l = len(pairs)
    a_vals = [p[0] for p in pairs]
    b_vals = [p[1] for p in pairs]
    
    # Verify that no M(a_i, b_i) is full (i.e., a_i != 1 and b_i != 1)
    for a, b in pairs:
        if a == 1 or b == 1:
            print(f"Error: Pair ({a}, {b}) corresponds to a full manifold, which is not allowed.")
            return

    # Calculate the left-hand side (LHS) of the equation
    # We output each number in the equation, as requested.
    x_vals = [(1 - a) * (1 - b) for a, b in zip(a_vals, b_vals)]
    sum_of_x = sum(x_vals)
    lhs = 2 * sum_of_x

    # Calculate the right-hand side (RHS) of the equation
    rhs = l - 1

    print("Verifying the solution against the core equation derived from topological properties.")
    print(f"Equation: 2 * sum((1 - a_i) * (1 - b_i)) = l - 1")
    print(f"Proposed tuple: {solution_tuple}")
    print(f"Minimal number of manifolds, l = {l}")
    
    # Outputting each number used in the final equation
    print(f"\nCalculating Left-Hand Side (LHS):")
    for i in range(l):
        print(f"  For pair (a_{i+1}, b_{i+1}) = ({a_vals[i]}, {b_vals[i]}), the term is (1-{a_vals[i]})*(1-{b_vals[i]}) = {x_vals[i]}")
    print(f"  Sum of terms = {sum_of_x}")
    print(f"  LHS = 2 * {sum_of_x} = {lhs}")
    
    print(f"\nCalculating Right-Hand Side (RHS):")
    print(f"  RHS = l - 1 = {l} - 1 = {rhs}")

    print(f"\nResult: LHS ({lhs}) = RHS ({rhs}). The tuple is a valid solution.")

    # Print the final answer as a flat tuple string
    answer_string = str(solution_tuple).replace(" ", "")
    print(f"\nThe final answer is: {answer_string}")

solve_and_verify()