import sys

def rsa_multiply_wuxing():
    """
    This script simulates RSA multiplication on a Wuxing-like architecture.
    It reads two large numbers p and q, multiplies them, and prints the result.
    The algorithm is optimized for the Wuxing constraints.
    """
    
    # --- Wuxing Architectural Constraints ---
    # int type is 5 decimal digits
    INT_DIGITS = 5
    # The numeric base for our array representation
    BASE = 10**INT_DIGITS
    # Maximum digits for inputs p and q
    P_Q_MAX_DIGITS = 100
    # Length of the array to hold p or q
    P_Q_ARRAY_LEN = P_Q_MAX_DIGITS // INT_DIGITS
    # Length of the array to hold the result o
    RESULT_ARRAY_LEN = P_Q_ARRAY_LEN * 2

    def parse_to_array(num_str):
        """Converts a number string into an array of 5-digit integers."""
        # Pad with leading zeros to make the length a multiple of INT_DIGITS
        num_str = num_str.zfill(P_Q_MAX_DIGITS)
        arr = []
        # Process the string from right to left in 5-digit chunks
        for i in range(P_Q_ARRAY_LEN):
            start = len(num_str) - (i + 1) * INT_DIGITS
            end = len(num_str) - i * INT_DIGITS
            chunk = num_str[start:end]
            arr.append(int(chunk))
        return arr

    try:
        p_str = sys.stdin.readline().strip()
        q_str = sys.stdin.readline().strip()
        
        if not p_str or not q_str:
            print("0")
            return

        # Convert input strings to our array representation
        p_arr = parse_to_array(p_str)
        q_arr = parse_to_array(q_str)
        
        # Initialize result array with zeros
        o_arr = [0] * RESULT_ARRAY_LEN

        # Perform multiplication using the schoolbook algorithm
        # This mirrors the C implementation for Wuxing
        for i in range(P_Q_ARRAY_LEN):
            carry = 0
            for j in range(P_Q_ARRAY_LEN):
                # The product of two 5-digit ints fits in a 10-digit long.
                # Python's integers handle this automatically.
                product = q_arr[i] * p_arr[j] + o_arr[i + j] + carry
                
                # Store the lower 5 digits of the result
                o_arr[i + j] = product % BASE
                
                # The upper digits are the new carry
                carry = product // BASE
                
            # Propagate the final carry from the inner loop
            if carry > 0:
                o_arr[i + P_Q_ARRAY_LEN] += carry
        
        # Format the result array back into a single number string for printing
        first_digit_idx = -1
        for i in range(RESULT_ARRAY_LEN - 1, -1, -1):
            if o_arr[i] != 0:
                first_digit_idx = i
                break
        
        if first_digit_idx == -1:
            print("0")
            return

        # Print the most significant chunk first
        result_str = str(o_arr[first_digit_idx])
        
        # Print the remaining chunks, padded with leading zeros
        for i in range(first_digit_idx - 1, -1, -1):
            result_str += f"{o_arr[i]:0{INT_DIGITS}d}"

        print(result_str)

    except Exception:
        # Fallback for any unexpected errors
        print("0")

rsa_multiply_wuxing()