import math

def decide_language_L(input_string):
    """
    This function acts as a decider for the language L.
    L contains all strings 'w' where the length of 'w' is a perfect square.
    It simulates the logic needed to determine if the Turing machine T would halt.
    """
    n = len(input_string)
    
    # Calculate the integer part of the square root of n
    k = int(math.sqrt(n))

    # Check if k squared equals n. This is a reliable way to check for a perfect square.
    is_perfect_square = (k * k == n)

    if is_perfect_square:
        print(f"Input string: '{input_string}'")
        print(f"Length of the string is n = {n}.")
        # The equation shows why n is a perfect square
        print(f"The Turing machine T halts because n is a perfect square: {k} * {k} = {n}.")
        print(f"Therefore, the string is in the decidable language L.")
    else:
        print(f"Input string: '{input_string}'")
        print(f"Length of the string is n = {n}.")
        print(f"The Turing machine T does not halt because n is not a perfect square.")
        print(f"Therefore, the string is not in the decidable language L.")
    
    print("-" * 40)

def main():
    """
    Main function to run the decider on a set of test strings.
    """
    test_strings = [
        "",           # length 0
        "a",          # length 1
        "abc",        # length 3
        "abcd",       # length 4
        "computer",   # length 8
        "recursion",  # length 9
        "decidability" # length 12
    ]

    print("Running the decider for language L on several test strings...\n")
    for s in test_strings:
        decide_language_L(s)

if __name__ == "__main__":
    main()