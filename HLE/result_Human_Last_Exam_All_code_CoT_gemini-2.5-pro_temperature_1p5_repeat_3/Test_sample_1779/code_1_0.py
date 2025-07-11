import math

def decide_language_of_perfect_squares(input_string):
    """
    Decides if a string w is in the language L = {w | |w| is a perfect square}.
    This function simulates a Turing machine that decides the language.
    It is guaranteed to halt for any input.
    """
    n = len(input_string)
    print(f"--- Analyzing input string: '{input_string}' ---")
    print(f"The length of the string is n = {n}.")

    if n < 0: # Should not happen for string length, but good practice
        print("Length is negative, which is not a perfect square. String is not in L.")
        return

    i = 0
    while True:
        square = i * i
        # The "final equation" step is to show the comparison being made.
        print(f"Checking i = {i}: is {i} * {i} = {square} equal to {n}?")

        if square == n:
            print(f"Yes! The final equation is {i} * {i} = {n}.")
            print(f"The length {n} is a perfect square. The string IS in the language L.")
            break
        elif square > n:
            print(f"No, {square} is greater than {n}. We have passed our number.")
            print(f"The length {n} is not a perfect square. The string is NOT in the language L.")
            break
        else: # square < n
            i += 1
            # Continue the loop

# Main part of the script to run examples
if __name__ == "__main__":
    # Example 1: A string with length 9 (a perfect square)
    w1 = "aaaaaaaaa"
    decide_language_of_perfect_squares(w1)

    print("\n" + "="*40 + "\n")

    # Example 2: A string with length 10 (not a perfect square)
    w2 = "bbbbbbbbbb"
    decide_language_of_perfect_squares(w2)

    print("\n" + "="*40 + "\n")

    # Example 3: A string with length 0 (a perfect square)
    w3 = ""
    decide_language_of_perfect_squares(w3)