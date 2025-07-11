import sys

def multiply_large_numbers():
    """
    This script simulates a C program on the Wuxing architecture to multiply
    two large integers, p and q.

    The Wuxing architecture has a 5-digit (5D) 'int' type. We leverage this
    by treating the input numbers as numbers in base 100,000. Each 'digit'
    in this base is a 5-digit number that fits into a Wuxing 'int'.
    """

    # Wuxing architecture constants
    CHUNK_SIZE = 5
    BASE = 10**CHUNK_SIZE

    def to_chunks(num_str):
        """Converts a number string into a list of base-100000 'digits'."""
        # Pad with leading zeros to make length a multiple of CHUNK_SIZE
        if len(num_str) % CHUNK_SIZE != 0:
            pad_len = CHUNK_SIZE - (len(num_str) % CHUNK_SIZE)
            num_str = '0' * pad_len + num_str
        
        # Handle empty or zero string
        if not num_str or int(num_str) == 0:
            return [0]

        chunks = []
        for i in range(0, len(num_str), CHUNK_SIZE):
            chunks.append(int(num_str[i:i+CHUNK_SIZE]))
        
        # Reverse so least significant chunk is at index 0
        return chunks[::-1]

    def format_result(res_chunks):
        """Converts a list of chunks back to a number string."""
        # Find the most significant non-zero chunk
        first_digit_idx = len(res_chunks) - 1
        while first_digit_idx > 0 and res_chunks[first_digit_idx] == 0:
            first_digit_idx -= 1
        
        # Format the first chunk without padding
        result_str = str(res_chunks[first_digit_idx])
        
        # Format the rest of the chunks with leading zero padding
        for i in range(first_digit_idx - 1, -1, -1):
            result_str += f"{res_chunks[i]:0{CHUNK_SIZE}d}"
            
        return result_str

    # --- Main Logic ---
    try:
        # Read inputs p and q
        p_str = sys.stdin.readline().strip()
        q_str = sys.stdin.readline().strip()

        if not p_str or not q_str:
            raise ValueError("Input cannot be empty.")

    except (IOError, ValueError) as e:
        print(f"Error reading input: {e}", file=sys.stderr)
        return

    # Convert strings to chunks (simulating arrays of 'int' on Wuxing)
    p_chunks = to_chunks(p_str)
    q_chunks = to_chunks(q_str)
    
    len_p = len(p_chunks)
    len_q = len(q_chunks)
    
    # Initialize result array (simulating 'int o[40]' on Wuxing)
    # The result of p * q has at most len(p) + len(q) digits.
    # In our chunked representation, it has at most len_p + len_q chunks.
    res_chunks = [0] * (len_p + len_q)

    # Schoolbook multiplication using chunks
    for i in range(len_q):
        if q_chunks[i] == 0:
            continue
        for j in range(len_p):
            # C equivalent: long prod = (long)p_chunks[j] * q_chunks[i];
            prod = p_chunks[j] * q_chunks[i]
            
            # Add product to the result array at the correct position and handle carries
            k = i + j
            
            # The sum can be large, but the carry propagation handles it
            current_val = res_chunks[k] + prod
            res_chunks[k] = current_val % BASE
            carry = current_val // BASE
            
            # Propagate the carry
            k += 1
            while carry > 0:
                # Ensure we don't go out of bounds (should not happen with correct sizing)
                if k >= len(res_chunks):
                    res_chunks.append(0)
                
                current_val = res_chunks[k] + carry
                res_chunks[k] = current_val % BASE
                carry = current_val // BASE
                k += 1

    # Format the result for printing
    result_str = format_result(res_chunks)
    
    # Print the final equation as requested
    print(f"{p_str} * {q_str} = {result_str}")


if __name__ == "__main__":
    multiply_large_numbers()