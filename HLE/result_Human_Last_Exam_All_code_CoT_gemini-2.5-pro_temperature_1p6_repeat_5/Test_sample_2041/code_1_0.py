def solve_shallow_functions():
    """
    Calculates the number of extensionally distinct functions induced by "shallow" expressions.

    The plan is as follows:
    1.  Analyze functions that can be formed WITHOUT using the variable p (p-independent).
    2.  Analyze functions that can be formed WITH the variable p, subject to the shallow constraint.
    3.  Sum the results, as the two sets of functions are disjoint.
    """

    print("### Step 1: Counting `p`-independent functions ###")
    print("A function F(p, x) is p-independent if its body `e` is a Bool expression built only from x:X.")
    print("A Bool is encoded as `λt:X. λf:X. body`, where `body` has type X.")
    print("To construct `e`, the available variables of type X for the `body` are `t`, `f`, and the free variable `x`.")
    
    # The choices for the body are the variables t, f, and x.
    # 1. body = t  --> λt.λf.t (True)
    # 2. body = f  --> λt.λf.f (False)
    # 3. body = x  --> λt.λf.x (A Bool that always resolves to the value of x)
    num_p_independent_functions = 3
    
    print(f"Number of choices for the body of `e`: {num_p_independent_functions} (can be t, f, or x).")
    print(f"This yields {num_p_independent_functions} distinct functions that ignore `p`.\n")

    print("### Step 2: Counting `p`-dependent functions ###")
    print("A shallow p-dependent function F(p, x) must have a body `e` of the form `p(arg1)(arg2)(arg3)`.")
    print("The arguments `arg1`, `arg2`, `arg3` must be of type PX (X -> Bool) and cannot depend on `p`.")

    print("\nFirst, we count the number of valid `p`-independent arguments (predicates of type PX).")
    print("An argument `arg` has the form `λy:X. b`, where `b` is a Bool built from `x:X` and `y:X`.")
    print("The Bool `b` is `λt:X. λf:X. body`, where `body` has type X.")
    print("Available variables of type X for this body are `t`, `f`, `y` (from λy), and `x` (global).")

    # The choices for the body of the argument's Bool are t, f, y, and x.
    # 1. body = t --> λy.λt.λf.t (Predicate that always returns True)
    # 2. body = f --> λy.λt.λf.f (Predicate that always returns False)
    # 3. body = y --> λy.λt.λf.y (Predicate that returns a Bool resolving to its input y)
    # 4. body = x --> λy.λt.λf.x (Predicate that returns a Bool resolving to the global x)
    num_p_independent_predicates = 4
    
    print(f"Number of choices for the body of the argument: {num_p_independent_predicates} (can be t, f, y, or x).")
    print(f"This yields {num_p_independent_predicates} distinct p-independent predicates to use as arguments for `p`.\n")

    print("Now, we count the number of ways to form `p(arg1)(arg2)(arg3)`.")
    print("`p` takes 3 arguments, and for each argument, we can choose any of the `p`-independent predicates.")
    
    num_p_arguments = 3
    num_p_dependent_functions = num_p_independent_predicates ** num_p_arguments
    
    print(f"Number of distinct functions = (choices for arg1) * (choices for arg2) * (choices for arg3)")
    print(f"Calculation: {num_p_independent_predicates}^{num_p_arguments} = {num_p_dependent_functions}")
    print(f"This yields {num_p_dependent_functions} distinct functions that depend on `p`.\n")

    print("### Step 3: Total Count ###")
    print("The total number of shallow functions is the sum of the p-independent and p-dependent functions.")
    
    total_functions = num_p_independent_functions + num_p_dependent_functions
    
    print(f"Total = (p-independent count) + (p-dependent count)")
    print(f"Total = {num_p_independent_functions} + {num_p_dependent_functions} = {total_functions}")

    print(f"\n<<<{total_functions}>>>")

solve_shallow_functions()