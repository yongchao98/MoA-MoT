import math

def check_string_length(length_w):
    """
    This function simulates a decider for the language L.
    It takes the length of a string 'w' and determines if it's a perfect square.
    The underlying algorithm is guaranteed to halt for any non-negative integer input.

    Args:
        length_w: A non-negative integer representing the length of string 'w'.
    """
    if not isinstance(length_w, int) or length_w < 0:
        print(f"Error: The length must be a non-negative integer. Received: {length_w}")
        return

    # A simple and efficient way to check for a perfect square.
    # We find the integer part of the square root.
    root = int(math.sqrt(length_w))

    # Then, we check if squaring this root gives back the original number.
    if root * root == length_w:
        # If it does, the length is a perfect square. The string is in the language L.
        print(f"ACCEPT: The length {length_w} is a perfect square.")
        # Printing each number in the final equation as requested.
        print(f"{length_w} = {root} * {root}")
    else:
        # If not, the length is not a perfect square. The string is not in L.
        print(f"REJECT: The length {length_w} is not a perfect square.")

# --- Demonstration of the Decider's Logic ---
# Test Case 1: A string of length 49 (which is in L)
print("Testing with length 49...")
check_string_length(49)
print("-" * 20)

# Test Case 2: A string of length 50 (which is not in L)
print("Testing with length 50...")
check_string_length(50)
print("-" * 20)

# Test Case 3: A string of length 0 (which is in L)
print("Testing with length 0...")
check_string_length(0)
print("-" * 20)

# Test Case 4: A string of length 1 (which is in L)
print("Testing with length 1...")
check_string_length(1)
print("-" * 20)
