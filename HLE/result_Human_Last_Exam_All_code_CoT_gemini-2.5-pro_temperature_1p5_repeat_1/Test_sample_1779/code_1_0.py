import math

def simulate_decider_for_perfect_square_length(input_string):
    """
    This function simulates a decider for the language L where a string w
    is in L if and only if its length is a perfect square.

    It prints the step-by-step process of the decision.
    """
    w = input_string
    n = len(w)
    
    print(f"--- Analyzing string: \"{w}\" ---")
    print(f"Length of the string is n = {n}")

    # A decider must halt on all inputs. We'll simulate this with a loop
    # that is guaranteed to terminate.
    k = 0
    while True:
        square = k * k
        
        # In each step, we check the condition
        print(f"Checking k={k}: is {k} * {k} == {n}? Current square is {square}.")

        if square == n:
            # We found the integer square root. The length is a perfect square.
            # The string is in the language.
            print(f"Result: ACCEPT. The length is a perfect square.")
            print(f"Final Equation: {k} * {k} = {n}")
            break
        elif square > n:
            # The square has exceeded the length n. It's not a perfect square.
            # The string is not in the language.
            print(f"Result: REJECT. The length is not a perfect square.")
            break
        
        # If we haven't accepted or rejected, we continue the process.
        k += 1
    print("-" * 30)

# --- Simulation with example strings ---

# Example 1: A string whose length (9) is a perfect square
string1 = "It works!"
simulate_decider_for_perfect_square_length(string1)

# Example 2: A string whose length (12) is not a perfect square
string2 = "Hello World!"
simulate_decider_for_perfect_square_length(string2)

# Example 3: The empty string, whose length (0) is a perfect square
string3 = ""
simulate_decider_for_perfect_square_length(string3)
