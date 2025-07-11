import math

def demonstrate_decidability(input_string):
    """
    This function demonstrates the algorithm to decide the language L.
    L = {w | len(w) is a perfect square}.
    The existence of this always-halting algorithm proves L is decidable.
    """
    n = len(input_string)
    print(f"Analyzing input string w of length n = {n}")
    print("-" * 30)
    print("Algorithm to decide if n is a perfect square:")

    k = 0
    while True:
        k_squared = k * k
        print(f"Checking k={k}: {k} * {k} = {k_squared}")

        if k_squared == n:
            print(f"\nResult: Found a match!")
            print(f"The final equation is: {k} * {k} = {n}")
            print(f"The length of the string ({n}) is a perfect square.")
            print("Therefore, the string w is in the language L.")
            break
        elif k_squared > n:
            print(f"\nResult: No match found.")
            print(f"The deciding equation is: {k} * {k} = {k_squared}, which is greater than n={n}.")
            print(f"The length of the string ({n}) is not a perfect square.")
            print("Therefore, the string w is NOT in the language L.")
            break
        k += 1

    print("-" * 30)
    print("Conclusion: Since this algorithm is guaranteed to halt for any input string,")
    print("the language L is decidable.")
    print("So, a decidable language L as described in the question exists.")

# Example 1: A string with a length that is a perfect square.
demonstrate_decidability("a" * 16)

print("\n" + "="*50 + "\n")

# Example 2: A string with a length that is not a perfect square.
demonstrate_decidability("a" * 18)
