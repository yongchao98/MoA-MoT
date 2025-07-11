import math

def demonstrates_decider_for_L(input_string):
    """
    This function simulates a decider for the language L, where L is the set of
    strings whose length is a perfect square. The decider always halts.
    """
    n = len(input_string)
    
    print(f"The Turing Machine T receives the input string w = '{input_string}'.")
    print(f"The language L contains all strings w where T halts.")
    print(f"T halts if and only if the length of w is a perfect square.")
    print(f"Step 1: Calculate the length of w. |w| = {n}.")
    print(f"Step 2: Check if {n} is a perfect square.")

    # A length cannot be negative, but handle the case n=0 explicitly.
    if n == 0:
        is_perfect_sq = True
        root = 0
    else:
        # Calculate the integer part of the square root
        root = int(math.sqrt(n))
        # Check if squaring the integer root gives back the original number
        if root * root == n:
            is_perfect_sq = True
        else:
            is_perfect_sq = False

    if is_perfect_sq:
        print(f"Result: The length {n} is a perfect square.")
        print(f"Equation: {root} * {root} = {n}")
        print("This means T halts on w, so w is in L.")
    else:
        print(f"Result: The length {n} is not a perfect square.")
        print("This means T does not halt on w, so w is not in L.")
    
    print("Since this check can be performed by an algorithm that always halts, the language L is decidable.")

# Let's run the demonstration with an example string of length 16.
demonstrates_decider_for_L("a_perfect_square")
