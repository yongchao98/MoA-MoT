import sys

def large_multiply(p_str, q_str):
    """
    Multiplies two large numbers represented as strings using the
    schoolbook algorithm, simulating the Wuxing architecture's
    decimal-based digit-by-digit computation.
    """
    # On Wuxing, p, q, o would be char arrays. In Python, we use lists of ints.
    # Convert string digits to ints and reverse for easier math (LSB at index 0)
    # e.g., "123" -> [3, 2, 1]
    if p_str == '0' or q_str == '0':
        return '0'
        
    p_digits = [int(d) for d in p_str][::-1]
    q_digits = [int(d) for d in q_str][::-1]
    
    len_p = len(p_digits)
    len_q = len(q_digits)
    
    # The result 'o' will have at most len_p + len_q digits.
    # Initialize with zeros. This corresponds to `char o[200] = {0};`
    o_digits = [0] * (len_p + len_q)
    
    # 1. Multiply digit by digit without carrying over yet.
    #    On Wuxing, each 'product' would be stored in a 5D 'int' register.
    for i in range(len_p):
        for j in range(len_q):
            o_digits[i+j] += p_digits[i] * q_digits[j]

    # 2. Normalize the result array by handling carries.
    #    This is a separate loop for clarity.
    carry = 0
    for i in range(len(o_digits)):
        # The sum fits in a Wuxing 'int' register
        total = o_digits[i] + carry
        o_digits[i] = total % 10  # The new digit
        carry = total // 10       # The carry to the next position
        
    # Remove leading zeros from the result
    while len(o_digits) > 1 and o_digits[-1] == 0:
        o_digits.pop()
        
    # Convert result list back to a string and print the equation
    o_str = "".join(map(str, o_digits[::-1]))
    
    return o_str

def main():
    """
    Main function to read inputs and print the full multiplication equation.
    """
    try:
        # Read the two large numbers from stdin.
        p_str = input().strip()
        q_str = input().strip()

        # Validate input to be non-negative integers
        if not p_str.isdigit() or not q_str.isdigit():
            print("Error: Inputs must be non-negative integers.", file=sys.stderr)
            return

        o_str = large_multiply(p_str, q_str)
        
        # As requested, print the final equation including all numbers.
        print(f"{p_str} * {q_str} = {o_str}")

    except (IOError, ValueError) as e:
        print(f"An error occurred: {e}", file=sys.stderr)

if __name__ == "__main__":
    main()