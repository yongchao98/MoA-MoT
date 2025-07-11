import math

def solve_turing_problem(input_string: str):
    """
    This function determines if a given string belongs to the language L,
    where L is the set of strings whose length is a perfect square.
    This simulates a "decider" for the language L.

    The existence of this decider, which always halts and gives a correct
    yes/no answer, proves that the language L is decidable.
    """
    
    # Step 1: The decider calculates the length of the input string w.
    length = len(input_string)
    print(f"The input string is '{input_string}'.")
    print(f"Step 1: Calculate the length of the string. Length = {length}.")

    # Step 2: The decider checks if the length is a perfect square.
    # This is a finite computation that always terminates.
    print(f"Step 2: Check if the length ({length}) is a perfect square.")
    
    is_perfect = False
    # A negative length is not possible for a string, but we handle it for completeness.
    if length >= 0:
        # Calculate the integer part of the square root.
        root = int(math.sqrt(length))
        # Check if the square of the root equals the original number.
        if root * root == length:
            is_perfect = True
            
    # Step 3: The decider halts and provides a definitive answer.
    print("Step 3: Halt and output the result.")
    if is_perfect:
        print(f"Result: YES. The length {length} is a perfect square.")
        # As requested, output the numbers in the final equation.
        print(f"Equation: {root} * {root} = {length}")
        print("The Turing machine T halts on this input, so the string is in language L.")
    else:
        print(f"Result: NO. The length {length} is not a perfect square.")
        print("The Turing machine T does not halt on this input, so the string is not in language L.")

    print("\nConclusion: A decidable language L exists because we just demonstrated an algorithm that always halts and correctly decides membership in L.")


# --- Main execution ---
# Example 1: A string whose length IS a perfect square.
print("--- Example 1 ---")
solve_turing_problem("123456789") # Length is 9

print("\n" + "="*50 + "\n")

# Example 2: A string whose length IS NOT a perfect square.
print("--- Example 2 ---")
solve_turing_problem("hello") # Length is 5
