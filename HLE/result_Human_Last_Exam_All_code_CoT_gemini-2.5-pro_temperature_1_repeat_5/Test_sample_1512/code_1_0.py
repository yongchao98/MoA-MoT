import sys

def solve():
    """
    This script solves the large integer multiplication problem based on the
    constraints of the hypothetical Wuxing computer architecture.
    """

    # --- Wuxing Architecture Simulation Constants ---
    # An 'int' is 5 decimal digits (0-99999). We use this as our base.
    BASE = 100000
    # The number of digits in one 'int' chunk.
    CHUNK_SIZE = 5

    # --- Helper Functions ---

    def string_to_int_array(num_str):
        """
        Converts a string of digits into a list of integer chunks, with the
        least significant chunk first (little-endian). This simulates storing
        a large number in an array of 'int's on Wuxing.
        """
        if not num_str.isdigit() or int(num_str) == 0:
            return [0]
        
        chunks = []
        # Process the string from right to left in chunks of CHUNK_SIZE
        for i in range(len(num_str), 0, -CHUNK_SIZE):
            start_index = max(0, i - CHUNK_SIZE)
            chunk = int(num_str[start_index:i])
            chunks.append(chunk)
        return chunks

    def int_array_to_string(arr):
        """
        Converts a list of integer chunks back to a full number string.
        This simulates printing the large integer result from the array.
        """
        # Find the most significant non-zero chunk to avoid leading zeros in the output
        try:
            first_chunk_index = next(i for i in range(len(arr) - 1, -1, -1) if arr[i] != 0)
        except StopIteration:
            return "0"

        # Format the result string. The first chunk is not padded.
        result_str = str(arr[first_chunk_index])
        
        # Subsequent chunks are padded with leading zeros to CHUNK_SIZE.
        for i in range(first_chunk_index - 1, -1, -1):
            result_str += f"{arr[i]:0{CHUNK_SIZE}d}"
            
        return result_str

    # --- Main Multiplication Logic ---

    # Read input strings for p and q
    p_str = sys.stdin.readline().strip()
    q_str = sys.stdin.readline().strip()

    # Convert strings to our chunked array representation
    p_arr = string_to_int_array(p_str)
    q_arr = string_to_int_array(q_str)

    p_len = len(p_arr)
    q_len = len(q_arr)
    
    # The result 'o' will have at most p_len + q_len chunks.
    # This corresponds to: int o[p_len + q_len] = {0};
    o_len = p_len + q_len
    o_arr = [0] * o_len

    # The core multiplication algorithm, simulating the Wuxing process.
    # This is an efficient schoolbook multiplication method.
    for j in range(q_len):
        if q_arr[j] == 0:
            continue
        
        # 'carry' would be stored in a Wuxing 'long' register (10D).
        carry = 0
        
        for i in range(p_len):
            # 'prod' would be calculated using a 'long' to hold the intermediate value.
            # prod = p[i] * q[j] + o[i+j] + carry
            prod = p_arr[i] * q_arr[j] + o_arr[i+j] + carry
            o_arr[i+j] = prod % BASE
            carry = prod // BASE
            
        # The final carry for this pass is placed in the next position.
        if j + p_len < o_len:
            o_arr[j + p_len] += carry

    # Convert the result array back to a string for printing
    o_str = int_array_to_string(o_arr)

    # Print the final equation as requested
    print(f"{p_str} * {q_str} = {o_str}")

solve()