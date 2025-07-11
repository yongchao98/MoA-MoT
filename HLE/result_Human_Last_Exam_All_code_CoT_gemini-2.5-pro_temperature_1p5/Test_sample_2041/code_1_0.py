def solve():
    """
    This function explains and calculates the number of distinct functions
    induced by shallow expressions 'e'.
    """

    # Let's represent the building blocks of our functions symbolically.
    # The variable 'p' represents a term of type PPPX.
    # The variable 'x' represents a term of type X.

    # According to the analysis, shallow expressions are formed from a limited set of operations.

    # Case 1: The expression 'e' is p-free.
    # It must be a closed term, so we have two possibilities:
    p_free_expressions = ["True", "False"]
    num_p_free = len(p_free_expressions)
    print(f"Functions from p-free expressions:")
    for i, expr in enumerate(p_free_expressions):
        print(f"  {i+1}. λp.λx. {expr}")
    print("-" * 20)

    # Case 2: The expression 'e' contains 'p'.
    # Under a strict interpretation of 'shallow', 'e' has the form p(q) or Not(p(q)),
    # where 'q' is a p-free term of type PPX built from 'x'.

    # There are 4 such 'q' terms we can construct:
    q_terms = [
        "λg. g(x)        (evaluates predicate g at x)",
        "λg. Not(g(x))   (evaluates g at x and negates)",
        "λg. True        (always returns True)",
        "λg. False       (always returns False)"
    ]
    num_q_terms = len(q_terms)
    print("Four p-free arguments 'q' for p:")
    for i, q_term in enumerate(q_terms):
        print(f"  q{i+1} = {q_term}")
    print("-" * 20)
    
    # For each 'q', we can form two functions: one applying p, one applying Not(p(...)).
    p_dependent_expressions = []
    for i in range(num_q_terms):
        p_dependent_expressions.append(f"p(q{i+1})")
        p_dependent_expressions.append(f"Not(p(q{i+1}))")

    num_p_dependent = len(p_dependent_expressions)
    print("Functions from p-dependent expressions:")
    base_idx = num_p_free
    for i, expr in enumerate(p_dependent_expressions):
         print(f"  {base_idx+i+1}. λp.λx. {expr}")
    print("-" * 20)


    # The total number of distinct functions is the sum from both cases.
    total_functions = num_p_free + num_p_dependent

    print(f"Total number of distinct functions = {num_p_free} (p-free) + {num_p_dependent} (p-dependent) = {total_functions}")

solve()