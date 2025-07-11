import math

def decider_for_L(w: str):
    """
    This function acts as a decider for the language L.
    L = {w | the length of w is a perfect square}.
    It simulates the decider algorithm and prints its steps.
    A decider must always halt, which this function does.
    """
    n = len(w)
    print(f"--- Running decider for input string w = '{w}' ---")
    print(f"Step 1: Calculate the length of w. |w| = {n}.")
    print(f"Step 2: Check if n = {n} is a perfect square.")

    if n < 0: # Length can't be negative, but for robustness
        print(f"Result: Length {n} is not a perfect square.")
        print(f"Therefore, '{w}' is NOT in language L.\n")
        return

    # We test k = 0, 1, 2, ...
    k = 0
    while True:
        square = k * k
        print(f"  - Trying k = {k}. Is {k} * {k} = {square} equal to {n}?")

        if square == n:
            print(f"  - Yes! Found a match.")
            print(f"Result: Length {n} is a perfect square, as shown by the equation:")
            print(f"{n} = {k} * {k}")
            print(f"Therefore, '{w}' IS in language L.\n")
            return
        
        if square > n:
            print(f"  - No, and {square} > {n}. No larger k will work.")
            print(f"Result: Length {n} is not a perfect square.")
            print(f"Therefore, '{w}' is NOT in language L.\n")
            return
        
        k += 1

# Let's test the decider with some examples.
# A string of length 9 (which is 3*3), so it IS in L.
decider_for_L("perfectsq") 
# A string of length 7 (not a perfect square), so it is NOT in L.
decider_for_L("helloai")
# A string of length 0 (which is 0*0), so it IS in L.
decider_for_L("")
# A string of length 1 (which is 1*1), so it IS in L.
decider_for_L("a")
