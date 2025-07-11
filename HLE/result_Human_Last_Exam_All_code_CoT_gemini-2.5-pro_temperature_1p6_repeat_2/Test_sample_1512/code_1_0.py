import sys

def rsa_multiply():
    """
    This function reads two large numbers as strings, multiplies them,
    and prints the result. This demonstrates the core multiplication task
    described for the Wuxing architecture.
    
    The problem asks us to conceptualize an optimized C program for the
    Wuxing architecture and calculate the memory usage for the input
    and output numbers.

    - p: max 100 digits, stored as int p[20] (20 * 5D = 100D)
    - q: max 100 digits, stored as int q[20] (20 * 5D = 100D)
    - o: max 200 digits, stored as int o[40] (40 * 5D = 200D)

    The minimized total memory for p, q, and o is 100 + 100 + 200 = 400D.
    """
    try:
        # Read the first large number p from standard input
        p_str = sys.stdin.readline().strip()
        if not p_str:
            return

        # Read the second large number q from standard input
        q_str = sys.stdin.readline().strip()
        if not q_str:
            return

        # Convert the string representations to integers. Python handles
        # arbitrarily large integers automatically.
        p_val = int(p_str)
        q_val = int(q_str)
        
        # Perform the multiplication p * q to get o
        o_val = p_val * q_val
        
        # The prompt asks to "output each number in the final equation"
        # and provides an output format that is just "o". We will follow
        # the specified output format by printing only the result 'o'.
        print(o_val)

    except (ValueError, IndexError):
        # Handle cases with invalid input
        print("Please provide two valid integers, each on a new line.")

if __name__ == "__main__":
    rsa_multiply()
