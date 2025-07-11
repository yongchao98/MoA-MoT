import math

def check_string_for_perfect_square_length(w: str):
    """
    This function acts as a decider for the language L = {w | |w| is a perfect square}.
    It checks if the length of the input string 'w' is a perfect square and prints the result.
    """
    n = len(w)
    
    # A simple way to check for a perfect square is to take the square root,
    # round it to the nearest integer, and see if that integer squared equals the original number.
    
    if n < 0:
        # This case is theoretical as string length is non-negative.
        print(f'Input string: "{w}"\nLength: {n}')
        print("Error: Length is negative. Not in the language.\n")
        return

    if n == 0:
        print(f'Input string: "{w}"\nLength: {n}')
        # Show the equation
        print(f"Checking: 0 * 0 = 0")
        print(f"Result: The length {n} is a perfect square. The string is in the language L.\n")
        return
        
    # Calculate the integer part of the square root
    k = int(math.sqrt(n))
    
    print(f'Input string: "{w}"\nLength: {n}')

    # Check if k*k equals n
    if k * k == n:
        # Show the final equation
        print(f"Checking: {k} * {k} = {k*k}")
        print(f"Result: The length {n} is a perfect square. The string is in the language L.\n")
    else:
        # Show the equations that bound n
        print(f"Checking: {k} * {k} = {k*k}")
        next_k = k + 1
        print(f"Checking: {next_k} * {next_k} = {next_k*next_k}")
        print(f"Result: The length {n} is not a perfect square. The string is NOT in the language L.\n")


# --- Main execution ---
# Test cases demonstrating the decider's logic
test_string_1 = "aaaaaaaaa"  # Length 9, a perfect square
test_string_2 = "hello"      # Length 5, not a perfect square
test_string_3 = ""           # Length 0, a perfect square
test_string_4 = "decidable"  # Length 9, a perfect square
test_string_5 = "TuringMachine" # Length 13, not a perfect square

check_string_for_perfect_square_length(test_string_1)
check_string_for_perfect_square_length(test_string_2)
check_string_for_perfect_square_length(test_string_3)
check_string_for_perfect_square_length(test_string_4)
check_string_for_perfect_square_length(test_string_5)