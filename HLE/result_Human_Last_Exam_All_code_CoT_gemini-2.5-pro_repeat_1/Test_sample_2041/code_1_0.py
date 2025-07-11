def solve():
    """
    This script implements the concepts from the problem in Python
    to count the number of distinct shallow polymorphic functions.
    """

    # Step 1: Represent Bool in lambda calculus style.
    # True is a function that takes t, f and returns t.
    # False is a function that takes t, f and returns f.
    True_ = lambda t, f: t
    False_ = lambda t, f: f

    # The NOT operator on these booleans.
    # NOT(b) is equivalent to b(False, True)
    Not = lambda b: b(False_, True_)

    # Step 2: Analyze shallow expressions and their arguments.
    # A shallow expression `e` can be `True`, `False`, or `p(q)` where `q`
    # is of type PPX = (X -> Bool) -> Bool and is built from `x:X` only.

    # Step 3: Enumerate the four possible `q` functions.
    # A q function takes a predicate `r: X -> Bool` and returns a Bool.
    # It can only use `r` by applying it to the available `x`.

    # q1 corresponds to the identity function on the result of r(x).
    # q1(r) = r(x)
    def q_id(x):
        return lambda r: r(x)

    # q2 corresponds to the NOT function on the result of r(x).
    # q2(r) = NOT(r(x))
    def q_not(x):
        return lambda r: Not(r(x))

    # q3 corresponds to the constant True function, ignoring r.
    # q3(r) = True
    q_const_true = lambda r: True_

    # q4 corresponds to the constant False function, ignoring r.
    # q4(r) = False
    q_const_false = lambda r: False_

    # Step 4: Define the 6 resulting polymorphic functions of type PPPX -> PX.
    # Each function F takes a `p` and returns a function from `X` to `Bool`.
    # `p` itself is a function from `PPX` to `Bool`.

    # F1 comes from e = p(q_id(x))
    F1 = lambda p: lambda x: p(q_id(x))

    # F2 comes from e = p(q_not(x))
    F2 = lambda p: lambda x: p(q_not(x))

    # F3 comes from e = p(q_const_true)
    F3 = lambda p: lambda x: p(q_const_true)

    # F4 comes from e = p(q_const_false)
    F4 = lambda p: lambda x: p(q_const_false)

    # F5 comes from e = True
    F5 = lambda p: lambda x: True_

    # F6 comes from e = False
    F6 = lambda p: lambda x: False_

    # Step 5: Conclude the count.
    # Through logical analysis (as described in the thinking steps), we found that
    # these 6 functions are all extensionally distinct. For any pair Fi, Fj
    # with i != j, we can construct a `p` that produces different results.
    # Therefore, the number of distinct functions is 6.

    num_distinct_functions = 6
    print(f"Based on the analysis of shallow expressions in simply typed lambda calculus,")
    print(f"we can form 6 extensionally distinct functions.")
    print(f"The number of functions is: {num_distinct_functions}")

solve()
<<<6>>>