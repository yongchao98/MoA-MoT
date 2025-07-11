import math

def decide_language_L(input_string):
    """
    Decides if a string 'w' belongs to the language L, where L = {w | |w| is a perfect square}.
    This function implements the decision algorithm.
    """
    n = len(input_string)
    print(f"--- Checking string \"{input_string}\" of length {n} ---")

    if n < 0: # Length cannot be negative
        is_member = False
    elif n == 0:
        is_member = True
        # Fulfilling the requirement to print the equation
        print("0 * 0 = 0")
    else:
        is_member = False
        # Iterate from 1 up to floor(sqrt(n)) + 1
        # The loop is guaranteed to terminate
        i = 1
        while i * i <= n:
            if i * i == n:
                # Found the integer square root, so n is a perfect square
                is_member = True
                # Fulfilling the requirement to print the equation
                print(f"{i} * {i} = {n}")
                break # Exit the loop
            i += 1

    if is_member:
        print(f"Result: The string is IN language L because its length {n} is a perfect square.")
    else:
        print(f"Result: The string is NOT in language L because its length {n} is not a perfect square.")
    print("-" * 30)


# --- Main execution ---
# Demonstrate the decider with a few example strings
# A string of length 9 (a perfect square)
w1 = "123456789"
decide_language_L(w1)

# A string of length 10 (not a perfect square)
w2 = "abcdefghij"
decide_language_L(w2)

# A string of length 16 (a perfect square)
w3 = "a" * 16
decide_language_L(w3)

# A string of length 0 (a perfect square)
w4 = ""
decide_language_L(w4)
