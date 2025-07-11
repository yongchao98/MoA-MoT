import math

def decider_for_L(w: str):
    """
    This function decides if a string 'w' belongs to the language L,
    where L = {w | len(w) is a perfect square}.
    It simulates the decision procedure.
    """
    n = len(w)
    print(f"Input string: '{w}'")
    print(f"Step 1: Calculate the length of the string, n = {n}")

    if n < 0:
        # Length cannot be negative, but for robustness:
        print("Length is negative, which is not a perfect square.")
        print("Result: NO (w is not in L)")
        return

    # Step 2: Find the integer part of the square root of n.
    k = math.isqrt(n) # Equivalent to int(n**0.5)
    print(f"Step 2: Find the integer square root of n, k = {k}")

    # Step 3: Check if k*k equals n. This is the 'final equation'.
    k_squared = k * k
    print(f"Step 3: Check the equation k * k = n")
    print(f"Calculation: {k} * {k} = {k_squared}")

    # Step 4: Make the decision.
    if k_squared == n:
        print(f"Result: YES. Since {k_squared} == {n}, the length is a perfect square, and w is in L.")
    else:
        print(f"Result: NO. Since {k_squared} != {n}, the length is not a perfect square, and w is not in L.")
    print("-" * 20)

# Main execution
print("Answering the question: Does there exist a decidable language L such that w in L iff T halts on w?")
print("The answer is YES, because we can construct an algorithm that always halts and decides L.\n")
print("Demonstrating the decision algorithm with examples:")
print("-" * 20)

# Example 1: A string whose length is a perfect square (9)
decider_for_L("perfectsq")

# Example 2: A string whose length is not a perfect square (10)
decider_for_L("notperfect")

# Example 3: A string of length 0 (0 is a perfect square, 0*0=0)
decider_for_L("")