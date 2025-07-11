import math

def decider_for_L(w: str):
    """
    This function decides if a string 'w' is in the language L, where L is the
    set of strings whose length is a perfect square.
    It simulates a decider Turing Machine by always halting and providing a
    'yes' or 'no' answer.
    """
    n = len(w)
    print(f"--- Checking string \"{w}\" ---")
    print(f"Input string is: '{w}'")
    print(f"Length of string |w| is: n = {n}")
    
    if n < 0:
        print("Length is negative, which is not possible. Rejecting.")
        return False

    k = 0
    while True:
        square = k * k
        print(f"Checking k = {k}: The equation is {k} * {k} = {square}")
        
        if square == n:
            print(f"Success! The length {n} is a perfect square.")
            print(f"Final Equation: {k} * {k} = {n}")
            print(f"Therefore, \"{w}\" is in the language L.")
            return True
        
        if square > n:
            print(f"Failed. {square} > {n}, so {n} cannot be a perfect square.")
            print(f"Therefore, \"{w}\" is NOT in the language L.")
            return False
            
        k += 1

# --- Main execution ---
# The question is: Does there exist a decidable language L such that w is in L
# if and only if T halts on w?
# The answer is Yes, because the language L itself is decidable.
# The code above is a demonstration of the decider for L.
# Let's run it with a few examples.

decider_for_L("#########")  # Length 9, is a perfect square
print("\n" + "="*30 + "\n")
decider_for_L("hello")     # Length 5, is not a perfect square
print("\n" + "="*30 + "\n")
decider_for_L("")           # Length 0, is a perfect square (0*0=0)
print("\n" + "="*30 + "\n")
decider_for_L("a")          # Length 1, is a perfect square (1*1=1)

print("\nAs shown, we can always determine membership in L in a finite number of steps.")
print("The procedure never loops infinitely. Therefore, the language L is decidable.")
print("The answer to the question 'Does there exist such a decidable language?' is Yes.")
