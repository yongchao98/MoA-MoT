import sys

# This script performs the multiplication of two large numbers.
# The memory calculation 'm' is derived from the constraints of the 
# hypothetical Wuxing architecture, as explained in the steps above.

# Read p and q from standard input.
# Assumes p and q are provided on two separate lines.
try:
    p_str = sys.stdin.readline().strip()
    q_str = sys.stdin.readline().strip()

    # Ensure inputs are not empty and are valid digits
    if not p_str.isdigit() or not q_str.isdigit():
        raise ValueError("Inputs must be non-negative integers.")

    # Convert strings to integers for calculation.
    # Python's int type handles arbitrarily large numbers.
    p_int = int(p_str)
    q_int = int(q_str)

    # Calculate the product.
    o_int = p_int * q_int

    # Convert the result back to a string for printing.
    o_str = str(o_int)

    # Per the instruction to "output each number in the final equation",
    # we print the full expression.
    print(f"{p_str} * {q_str} = {o_str}")

except (ValueError, IndexError) as e:
    print(f"Error: Invalid input. Please provide two integers on separate lines. Details: {e}", file=sys.stderr)
