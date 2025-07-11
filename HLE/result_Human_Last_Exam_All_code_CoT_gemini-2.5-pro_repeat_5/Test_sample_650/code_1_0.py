def imag_base_converter(n, c):
    """
    Implements the IMAG(n, c) algorithm to represent number n in base c.
    This function also prints the verification equation.

    Args:
        n (int): The non-negative integer to convert.
        c (int): The base to convert to (must be >= 2).

    Returns:
        None. Prints the result directly.
    """
    if not isinstance(n, int) or n < 0:
        print("Error: n must be a non-negative integer.")
        return
    if not isinstance(c, int) or c < 2:
        print("Error: c must be an integer >= 2.")
        return
        
    original_n = n
    
    # Step 1: Initialization
    i = 0
    x = n
    q = x // c
    a0 = x - q * c
    
    # Store digits in a list. The list will be built in reverse order.
    digits = [a0]
    
    # Step 2: Loop while q > 0
    while q > 0:
        # Step 2.1
        i += 1
        x = q
        q = x // c
        ai = x - q * c
        digits.append(ai)
    
    # The digits are in reverse order (a0, a1, a2, ...). Reverse it to get (ai, ..., a0).
    digits.reverse()
    
    # Step 3: Return the representation
    print(f"The base-{c} representation of {original_n} is: {tuple(digits)}")
    
    # Print the final verification equation
    equation_parts = []
    for power, digit in enumerate(reversed(digits)):
        equation_parts.append(f"{digit} * {c}^{power}")
    
    # Reverse the parts to match the order of the digits (highest power first)
    equation_parts.reverse()
    
    final_equation = f"{original_n} = " + " + ".join(equation_parts)
    print("Verification Equation:")
    print(final_equation)

# Example usage:
# Let's find the representation of n = 255 in base c = 16 (hexadecimal)
n_input = 255
c_input = 16
imag_base_converter(n_input, c_input)

print("-" * 20)

# Let's find the representation of n = 100 in base c = 2 (binary)
n_input = 100
c_input = 2
imag_base_converter(n_input, c_input)