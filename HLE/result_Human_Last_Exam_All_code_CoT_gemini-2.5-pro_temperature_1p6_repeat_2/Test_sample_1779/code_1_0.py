import math

def check_language_membership(input_string):
    """
    This function simulates a decider for the language L.
    The language L is defined as the set of strings w where T halts on w.
    T halts on w if and only if len(w) is a perfect square.

    The existence of this algorithm, which always halts and gives a correct
    yes/no answer, proves that the language L is decidable.
    """
    n = len(input_string)
    print(f"Input string: \"{input_string}\"")
    print(f"Length of the string is: {n}")

    is_square = False
    root = -1

    # An efficient algorithm to decide if n is a perfect square.
    # This check is a finite computation that always terminates.
    if n >= 0:
        root_candidate = int(math.sqrt(n))
        if root_candidate * root_candidate == n:
            is_square = True
            root = root_candidate

    # Output the result of the decision procedure.
    if is_square:
        print(f"Result: The length {n} is a perfect square.")
        # As requested, output the numbers in the final equation.
        print(f"Equation: {root} * {root} = {n}")
        print("Therefore, T halts on this string, and it BELONGS to language L.")
    else:
        print(f"Result: The length {n} is not a perfect square.")
        print("Therefore, T does not halt on this string, and it DOES NOT BELONG to language L.")

# --- Demonstration ---
print("Does there exist a decidable language L such that w in L iff T halts on w?")
print("Answer: Yes. The language L is the set of all strings whose length is a perfect square.")
print("The following code demonstrates the decider for this language.\n")

print("--- Test Case 1: A string in L ---")
check_language_membership("perfect!!")  # Length is 9, a perfect square

print("\n--- Test Case 2: A string not in L ---")
check_language_membership("not-perfect") # Length is 11, not a perfect square