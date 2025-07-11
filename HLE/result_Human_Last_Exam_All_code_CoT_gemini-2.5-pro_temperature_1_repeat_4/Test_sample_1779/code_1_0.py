import math

def decider_for_L(input_string):
    """
    This function decides if a string 'w' belongs to the language L,
    where L = {w | |w| is a perfect square}.
    It implements the decision algorithm described above.
    """
    n = len(input_string)
    print(f"Input string: '{input_string}'")
    print(f"Length of string (n): {n}")

    if n < 0:
        # Length cannot be negative, but as a safeguard.
        print(f"Length {n} is not a valid non-negative integer.")
        return

    # Let's find if there is an integer k such that k*k = n.
    k = 0
    while True:
        square = k * k
        print(f"Checking k = {k}, k*k = {square}")
        
        # This is the final equation check for acceptance
        if square == n:
            print(f"Success: The final equation is {k} * {k} = {n}.")
            print(f"The length {n} is a perfect square. The string is in language L.")
            break
        
        # This is the final condition check for rejection
        if square > n:
            print(f"Halt: {k} * {k} ({square}) > {n}.")
            print(f"The length {n} is not a perfect square. The string is NOT in language L.")
            break
            
        k += 1

# --- Test Cases ---
# Case 1: A string whose length is a perfect square (9)
decider_for_L("pythonisfun")
print("-" * 20)

# Case 2: A string whose length is not a perfect square (10)
decider_for_L("helloworld")
print("-" * 20)

# Case 3: A string of length 0 (0 is a perfect square, 0*0=0)
decider_for_L("")
print("-" * 20)

# Case 4: A string of length 1 (1 is a perfect square, 1*1=1)
decider_for_L("a")
