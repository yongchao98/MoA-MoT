def print_base_representation(n, c):
    """
    Implements the IMAG(n, c) algorithm to find and print the base-c representation of n.
    It prints each digit of the final representation as requested.

    Args:
        n (int): A non-negative integer.
        c (int): The base, an integer greater than or equal to 2.
    """
    # Validate the inputs
    if not isinstance(n, int) or not isinstance(c, int) or n < 0 or c < 2:
        print("Input error: Please provide integers n >= 0 and c >= 2.")
        return

    # Handle the special case where n is 0 as per the algorithm's logic.
    if n == 0:
        print(f"The base-{c} representation of 0 is: (0)")
        return

    # Step 1: Initialization
    i = 0
    x = n
    q = x // c
    a0 = x - q * c
    
    digits = [a0] # This will store the digits a_0, a_1, a_2, ...

    # Step 2: Loop while q > 0
    while q > 0:
        # Step 2.1
        i = i + 1
        x = q
        q = x // c
        ai = x - q * c
        digits.append(ai)
        # Step 2.2: od (end of loop iteration)

    # Step 3: Return the representation
    # The algorithm generates digits from least significant (a_0) to most significant (a_i).
    # For standard display, we reverse the list.
    digits.reverse()
    
    # Format the output to show the sequence of numbers (digits)
    # This fulfills the instruction "output each number in the final equation!".
    print(f"The base-{c} representation of {n} is:", end=" ")
    representation_str = ", ".join(map(str, digits))
    print(f"({representation_str})")

# --- Example Usage ---
# You can change these values to test the algorithm.
number_to_convert = 345
base = 16
print_base_representation(number_to_convert, base)

number_to_convert = 25
base = 2
print_base_representation(number_to_convert, base)