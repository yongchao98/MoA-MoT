def decider_for_perfect_square_language(string_length):
    """
    This function simulates a decider for the language L = {w | |w| is a perfect square}.
    It takes an integer representing the length of a string and determines
    if it's a perfect square. The algorithm is guaranteed to halt.
    
    Args:
        string_length (int): The length of the input string w. Must be non-negative.
    """
    if not isinstance(string_length, int) or string_length < 0:
        print("Error: Input length must be a non-negative integer.")
        return

    print(f"Running the decider for a string of length n = {string_length}...")
    
    i = 0
    while True:
        square = i * i
        # The prompt requires outputting each number in the final equation.
        # We will show each check that the algorithm performs.
        print(f"Checking counter i = {i}: The equation is {i} * {i} = {square}.")
        
        if square == string_length:
            print(f"Result: The length {string_length} IS a perfect square.")
            print("The algorithm halts and ACCEPTS the string.")
            return
        
        if square > string_length:
            print(f"Result: {square} is greater than {string_length}. The length {string_length} IS NOT a perfect square.")
            print("The algorithm halts and REJECTS the string.")
            return
            
        i += 1

# --- Main Execution ---
# Let's test with a length that is a perfect square.
input_length_1 = 25
decider_for_perfect_square_language(input_length_1)

print("\n" + "="*40 + "\n")

# Let's test with a length that is NOT a perfect square.
input_length_2 = 28
decider_for_perfect_square_language(input_length_2)
