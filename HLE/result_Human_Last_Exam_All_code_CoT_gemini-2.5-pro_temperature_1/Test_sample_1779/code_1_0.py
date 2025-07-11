import math

def simulate_decider_for_L(string_length):
    """
    This function simulates a Turing machine that decides the language L,
    where L is the set of strings whose length is a perfect square.

    A "decider" must halt on all inputs. This algorithm always halts.

    Args:
        string_length (int): The length of the input string 'w'.
    """
    print(f"--- Simulating decider for a string of length {string_length} ---")

    # Step 1: The length 'n' is given as string_length.
    n = string_length

    if n < 0:
        print(f"Error: Length cannot be negative. Rejecting.")
        return

    # Step 2: Determine if n is a perfect square.
    # We can find the integer closest to the square root of n.
    sqrt_n_int = int(math.sqrt(n))

    # Step 3: Check if squaring this integer gives us back n.
    # This is the "equation" check.
    square_of_sqrt = sqrt_n_int * sqrt_n_int
    
    print(f"Checking if length n = {n} is a perfect square.")
    print(f"The integer part of the square root of {n} is {sqrt_n_int}.")
    print(f"Final check equation: {sqrt_n_int} * {sqrt_n_int} = {square_of_sqrt}")

    # Step 4: Accept if it's a perfect square, reject otherwise.
    # The machine always halts after this finite check.
    if square_of_sqrt == n:
        print(f"Result: The equation holds. A string of length {n} is in L. Accepting.")
    else:
        print(f"Result: The equation does not hold. A string of length {n} is not in L. Rejecting.")

# --- Main execution ---
# Let's test two cases: one where the length is a perfect square and one where it is not.
length_in_L = 81
length_not_in_L = 50

simulate_decider_for_L(length_in_L)
print("\n")
simulate_decider_for_L(length_not_in_L)