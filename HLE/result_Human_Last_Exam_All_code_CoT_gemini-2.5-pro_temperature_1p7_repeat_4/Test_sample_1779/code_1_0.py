import math

def check_language_membership(input_string):
    """
    This function implements a decision algorithm for the language L.
    L = {w | |w| is a perfect square}.
    It prints the steps to determine if the input string's length is a perfect square.
    """
    n = len(input_string)
    print(f"The input string 'w' has a length of |w| = {n}.")
    print(f"We must decide if {n} is a perfect square.")
    print("-----------------------------------------")
    print("Algorithm steps:")

    # An efficient way to check for a perfect square.
    # We only need to check integers i up to sqrt(n).
    # If a perfect square root exists, it must be less than or equal to floor(sqrt(n)).
    if n < 0:
        is_member = False
        print(f"Length {n} is negative, so it cannot be a perfect square.")
    elif n == 0:
        is_member = True
        print(f"Checking i = 0: 0 * 0 = 0. This matches the length.")
    else:
        i = 1
        is_member = False
        while True:
            square = i * i
            print(f"Checking i = {i}: {i} * {i} = {square}")
            if square == n:
                print(f"Match found! {n} is a perfect square.")
                is_member = True
                break
            if square > n:
                print(f"Overshot! {square} > {n}, so {n} cannot be a perfect square.")
                is_member = False
                break
            i += 1

    print("-----------------------------------------")
    print("\nConclusion:")
    if is_member:
        print(f"The string 'w' of length {n} IS in the language L.")
    else:
        print(f"The string 'w' of length {n} IS NOT in the language L.")

    print("\nBecause we have an algorithm that always halts and correctly decides membership in L,")
    print("the language L is DECIDABLE.")
    print("Therefore, a decidable language L that satisfies the condition exists.")


# Let's test with a sample string of length 25 (a perfect square)
sample_string = "a" * 25
check_language_membership(sample_string)
