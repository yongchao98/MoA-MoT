import math

def check_if_length_is_perfect_square(w: str):
    """
    This function simulates a decider for the language L where L = {w | |w| is a perfect square}.
    It always halts and determines if the length of w is a perfect square.
    """
    n = len(w)
    print(f"Input string: '{w}'")
    print(f"Length of the string is n = {n}.")
    print("\nNow, we check if n is a perfect square by testing integers k = 0, 1, 2, ...")

    if n < 0:
        print(f"Length {n} is negative, so it cannot be a perfect square.")
        return False

    k = 0
    while True:
        square = k * k
        print(f"Checking: {k} * {k} = {square}")
        
        if square == n:
            print(f"\nSuccess! The length {n} is a perfect square because {k} * {k} = {n}.")
            print(f"Therefore, the Turing machine T halts on '{w}', and '{w}' is in the language L.")
            return True
        
        if square > n:
            print(f"\nFailure. The square {square} is greater than the length {n}.")
            print(f"No larger value of k will work, so {n} is not a perfect square.")
            print(f"Therefore, the Turing machine T does not halt on '{w}', and '{w}' is not in the language L.")
            return False
            
        k += 1

# --- Main Execution ---
# Example 1: A string whose length is a perfect square (16)
input_string_1 = "abcdefghijklmnop"
check_if_length_is_perfect_square(input_string_1)

print("\n" + "="*50 + "\n")

# Example 2: A string whose length is not a perfect square (10)
input_string_2 = "0123456789"
check_if_length_is_perfect_square(input_string_2)

print("\n" + "="*50 + "\n")
print("Conclusion: Yes, a decidable language L exists because we can create an algorithm (a decider) that always halts and correctly determines if a string's length is a perfect square.")
