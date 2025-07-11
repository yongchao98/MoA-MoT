def solve():
    """
    Calculates the number of extensionally distinct shallow polymorphic functions.

    The problem asks for the number of distinct functions F = λp.λx.e, where e is a "shallow"
    expression of type Bool.

    1.  A shallow expression e of type Bool = X->X->X in normal form is λt:X.λf:X.M, where M is a
        term of type X.

    2.  M's possible normal forms determine the possible functions e.
        The free variables available for M are {p, x, t, f}.

    3.  Case 1: M's head variable is x, t, or f.
        - M = x  => e = λt.λf.x
        - M = t  => e = true
        - M = f  => e = false
        This gives 3 functions that do not depend on p.
    """
    p_independent_functions = 3
    print(f"Number of functions not depending on p: {p_independent_functions}")
    print("These correspond to e being 'true', 'false', or 'λt.λf.x'.")

    """
    4.  Case 2: M's head variable is from an application of p.
        The shallow condition means p is only applied to p-free arguments A.
        M must be of the form p(A)(N1)(N2).

        - The argument A: There are 4 possible p-free arguments A that can be formed.
          A can be λq.q(x), λq.¬q(x), λq.true, or λq.false.
        
        - The arguments N1, N2: These must be terms of type X. The only available normal forms
          are x, t, and f. So, there are 3 choices for N1 and 3 choices for N2.
        
        - Total combinations: 4 choices for A * 3 choices for N1 * 3 choices for N2.
    """
    num_A_choices = 4
    num_N_choices = 3
    p_dependent_functions = num_A_choices * num_N_choices * num_N_choices
    print(f"\nNumber of functions depending on p:")
    print(f"Choices for argument A to p: {num_A_choices}")
    print(f"Choices for arguments N1, N2 to p(A): {num_N_choices}")
    print(f"Total p-dependent functions = {num_A_choices} * {num_N_choices} * {num_N_choices} = {p_dependent_functions}")

    """
    5.  Total distinct functions.
        The p-dependent functions are distinct from the p-independent ones. All functions within
        each group are also extensionally distinct.
    """
    total_functions = p_independent_functions + p_dependent_functions
    print(f"\nTotal number of distinct functions = {p_independent_functions} + {p_dependent_functions} = {total_functions}")

solve()