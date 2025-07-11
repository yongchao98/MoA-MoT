import sys

# This script demonstrates the multiplication of two large integers,
# as would be performed on the Wuxing architecture.
# It reads two numbers p and q from standard input and prints their product.

# Plan:
# 1. Read the first large integer, p, from standard input.
# 2. Read the second large integer, q, from standard input.
# 3. Use Python's built-in support for arbitrary-precision integers to
#    calculate the product, o = p * q.
# 4. Print the final equation, including the input numbers and the result,
#    as per the instruction: "output each number in the final equation".

try:
    # Read p and q from standard input, stripping any whitespace
    p_str = sys.stdin.readline().strip()
    q_str = sys.stdin.readline().strip()

    # Convert the string representations to integers. Python handles any size.
    p = int(p_str)
    q = int(q_str)

    # Calculate the product
    o = p * q

    # Print the full equation
    # The numbers in the equation are: p, q, and o.
    print(f"{p} * {q} = {o}")

except (ValueError, IndexError):
    # Handle cases with invalid or missing input
    print("Error: Please provide two valid integers on separate lines.", file=sys.stderr)
    sys.exit(1)