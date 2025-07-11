import math
import argparse

def check_if_length_is_in_language(length: int):
    """
    This function simulates a decider for the language L.
    L is the set of strings w where |w| is a perfect square.
    The function determines if a given length corresponds to a string in L.

    Args:
        length: An integer representing the length of the input string w.
    """
    if length < 0:
        print(f"Error: Length ({length}) cannot be negative.")
        return

    # A decider for L can be constructed. Here is its logic:
    # 1. It calculates the length of the input string, which is the number 'length'.
    #    This step always terminates.

    # 2. It checks if 'length' is a perfect square.
    #    This can be done by taking the integer part of the square root and squaring it.
    sqrt_val = int(math.sqrt(length))
    
    # The check itself is a finite, terminating operation.
    if sqrt_val * sqrt_val == length:
        # 3a. If it is a perfect square, the Turing Machine T halts.
        #     Therefore, the string is in the language L. The decider accepts.
        print(f"The length of the string is {length}.")
        print(f"Yes, a string of length {length} is in the language L.")
        print(f"This is because {length} is a perfect square.")
        # The final equation is printed as requested.
        print(f"Equation: {sqrt_val} * {sqrt_val} = {length}")
        print("Since a decider algorithm exists that always halts, the language L is decidable.")
    else:
        # 3b. If it's not a perfect square, T does not halt.
        #     The string is not in L. The decider rejects.
        print(f"The length of the string is {length}.")
        print(f"No, a string of length {length} is not in the language L.")
        print(f"This is because {length} is not a perfect square.")
        print("Since a decider algorithm exists that always halts, the language L is decidable.")

if __name__ == '__main__':
    # This allows the user to test with their own number from the command line.
    # Example usage: python your_script_name.py 16
    parser = argparse.ArgumentParser(
        description="Check if a string of a given length is in the language L where L = {w | |w| is a perfect square}."
    )
    parser.add_argument(
        'length', 
        type=int, 
        help='The length of the string to check.'
    )
    args = parser.parse_args()
    check_if_length_is_in_language(args.length)