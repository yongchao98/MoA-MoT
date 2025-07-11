import math

def check_string_membership_in_L():
    """
    This function simulates a decider for the language L.
    L is the set of strings w where |w| (length of w) is a perfect square.
    The decider determines if a string belongs to L based on its length.
    """
    try:
        # Step 1: Get the length 'n' of the input string 'w'.
        # We will ask the user for this length directly.
        input_str = input("Enter a non-negative integer for the length of the string w: ")
        n = int(input_str)

        if n < 0:
            print("Error: Length cannot be negative.")
            return

        # Step 2: Determine if n is a perfect square.
        # This is a terminating algorithm, guaranteeing the decider halts.
        
        # We find the integer part of the square root.
        root = math.isqrt(n)
        
        # Check if squaring the root gives us the original number.
        if root * root == n:
            # If it is a perfect square, the Turing Machine T halts,
            # and the string w is in the decidable language L.
            print(f"\nYES, a string of length {n} is in the language L.")
            print("The Turing machine T would halt on this string.")
            print("The deciding equation is:")
            print(f"{root} * {root} = {n}")
        else:
            # If it is not a perfect square, T does not halt,
            # and the string w is not in L.
            print(f"\nNO, a string of length {n} is not in the language L.")
            print("The Turing machine T would not halt on this string.")
            print(f"{n} is not a perfect square.")

    except ValueError:
        print("Invalid input. Please enter a valid non-negative integer.")

if __name__ == "__main__":
    check_string_membership_in_L()
