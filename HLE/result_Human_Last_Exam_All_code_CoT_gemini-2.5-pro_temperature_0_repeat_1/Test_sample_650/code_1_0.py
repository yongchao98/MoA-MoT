import math

def imag_representation(n, c):
    """
    Represents a non-negative integer n in base c using the IMAG(n, c) algorithm.
    It also prints the final representation and the corresponding polynomial equation.

    Args:
        n: An integer, where n >= 0.
        c: An integer, where c >= 2.
    """
    # Validate inputs
    if not isinstance(n, int) or n < 0:
        print("Error: n must be a non-negative integer.")
        return
    if not isinstance(c, int) or c < 2:
        print("Error: c must be an integer greater than or equal to 2.")
        return

    # Handle the special case of n = 0
    if n == 0:
        print(f"The base-{c} representation of 0 is: (0)")
        print("The final equation is: 0 = 0 * {c}^0")
        return

    # Step 1: Initialization
    x = n
    q = math.floor(x / c)
    a0 = x - q * c
    digits = [a0]

    # Step 2: Loop while q > 0
    while q > 0:
        # Step 2.1
        x = q
        q = math.floor(x / c)
        ai = x - q * c
        digits.append(ai)

    # The digits are generated from least significant to most significant.
    # Reverse the list for the standard representation (most significant first).
    digits.reverse()
    
    # Step 3: Output the results
    representation_str = "".join(map(str, digits))
    print(f"The base-{c} representation of {n} is: ({representation_str})")
    
    # Output the final equation showing each number
    equation_parts = []
    power = len(digits) - 1
    for digit in digits:
        equation_parts.append(f"{digit} * {c}^{power}")
        power -= 1
    
    print(f"The final equation is: {n} = {' + '.join(equation_parts)}")


# --- Example Usage ---
# Find the representation of n=25 in base c=2
# The expected result is (11001) because 25 = 1*16 + 1*8 + 0*4 + 0*2 + 1*1
print("--- Example 1: n=25, c=2 ---")
imag_representation(25, 2)

print("\n--- Example 2: n=255, c=16 ---")
# The expected result is (FF) because 255 = 15*16^1 + 15*16^0
# The algorithm will produce digits [15, 15]
imag_representation(255, 16)