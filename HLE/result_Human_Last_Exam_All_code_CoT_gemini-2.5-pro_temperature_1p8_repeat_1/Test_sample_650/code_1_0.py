def represent_in_base(n, c):
    """
    Implements the IMAG(n, c) algorithm to find the representation of
    a number n in base c. It then prints the resulting digits and the
    full equation to verify the result.

    Args:
        n (int): The number to convert, must be >= 0.
        c (int): The base, must be >= 2.
    """
    if not (isinstance(n, int) and n >= 0 and isinstance(c, int) and c >= 2):
        print("Error: Input must be integers n >= 0 and c >= 2.")
        return

    original_n = n

    # Handle the special case of n = 0
    if n == 0:
        print(f"The base-{c} representation of 0 is: (0)")
        # The prompt asks to output each number in the final equation.
        print(f"Final equation: 0 = 0 * {c}^0")
        return

    # Step 1: Initialization
    i = 0
    x = n
    q = x // c
    a0 = x - q * c
    digits = [a0]  # List to store the digits a_i, a_{i-1}, ...

    # Step 2: Loop while q > 0
    while q > 0:
        # Step 2.1
        i += 1
        x = q
        q = x // c
        ai = x - q * c
        digits.append(ai)
        # Step 2.2: od (end of loop body)

    # Step 3: The representation is (a_i, a_{i-1}, ..., a_0).
    # The digits list is currently [a_0, a_1, ... a_i], so we reverse it.
    digits.reverse()

    print(f"The base-{c} representation of {original_n} is: {tuple(digits)}")

    print("\nFinal equation, showing each number involved:")
    print(f"{original_n} = ", end="")
    
    num_digits = len(digits)
    for j, digit in enumerate(digits):
        power = num_digits - 1 - j
        # Print each part of the sum: digit * c^power
        print(f"{digit} * {c}^{power}", end="")
        if j < num_digits - 1:
            print(" + ", end="")
    print() # For a final newline

# Example usage with n = 2024, c = 8
n_example = 2024
c_example = 8
represent_in_base(n_example, c_example)