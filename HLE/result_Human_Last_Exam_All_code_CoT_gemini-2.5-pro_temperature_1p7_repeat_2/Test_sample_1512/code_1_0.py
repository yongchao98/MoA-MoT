import sys

def solve():
    """
    This function reads two large numbers as strings, multiplies them,
    and prints the result.
    The task is inspired by RSA encryption which requires multiplying two large prime numbers.
    While the prompt describes a hypothetical C program on a 'Wuxing' architecture,
    this Python script performs the actual computation requested.

    The details of the C program design and memory calculation are as follows:
    1. C Program Design:
       - Represent p (100-digit) and q (100-digit) as arrays of 20 5-digit integers.
       - Represent o (200-digit result) as an array of 40 5-digit integers.
       - Use schoolbook multiplication. An intermediate 10-digit 'long' would be
         used to hold the product of two 5-digit integers to avoid overflow.
    2. Memory Calculation (m):
       - Storage for p: 100 digits (100 D)
       - Storage for q: 100 digits (100 D)
       - Storage for o: 200 digits (200 D)
       - Total minimized memory (m) = 100 + 100 + 200 = 400 D.
    """
    try:
        # Read the two large numbers from standard input.
        p_str = sys.stdin.readline().strip()
        q_str = sys.stdin.readline().strip()

        # Convert the string representations to integers. Python handles arbitrarily large integers.
        p = int(p_str)
        q = int(q_str)

        # Perform the multiplication.
        o = p * q
        
        # The prompt instruction "output each number in the final equation!" is ambiguous.
        # To strictly follow the "Output: o" format, we print only the final product.
        print(o)

    except (IOError, ValueError) as e:
        print(f"An error occurred: {e}", file=sys.stderr)

solve()