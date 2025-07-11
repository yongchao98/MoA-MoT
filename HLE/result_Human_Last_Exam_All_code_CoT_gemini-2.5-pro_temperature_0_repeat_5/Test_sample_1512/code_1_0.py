import math

def solve_wuxing_multiplication():
    """
    This function simulates an optimized multiplication of two 100-digit numbers
    on the conceptual Wuxing computer architecture.
    """

    # Step 1: Define architecture and problem constants
    # In Wuxing, an 'int' is 5 decimal digits (5D). We use this as our base.
    CHUNK_DIGITS = 5
    BASE = 10**CHUNK_DIGITS

    # p and q are 100 digits each. The result 'o' can be up to 200 digits.
    P_DIGITS = 100
    Q_DIGITS = 100
    O_DIGITS = P_DIGITS + Q_DIGITS

    # Calculate the number of 'int' chunks needed for each number
    P_CHUNKS = P_DIGITS // CHUNK_DIGITS
    Q_CHUNKS = Q_DIGITS // CHUNK_DIGITS
    O_CHUNKS = O_DIGITS // CHUNK_DIGITS

    # Step 2: Helper functions to handle data representation
    def convert_str_to_chunks(s_num, num_digits, chunk_size):
        """
        Converts a string number into a list of integer chunks.
        The list is little-endian (least significant chunk at index 0),
        which simplifies the multiplication algorithm.
        e.g., "1234567890" -> [67890, 12345]
        """
        s_num = s_num.zfill(num_digits)
        chunks = []
        # Iterate from the end of the string (least significant part)
        for i in range(len(s_num), 0, -chunk_size):
            start = max(0, i - chunk_size)
            chunk_str = s_num[start:i]
            chunks.append(int(chunk_str))
        return chunks

    def convert_chunks_to_str(chunks, chunk_size):
        """
        Converts a little-endian list of chunks back to a string.
        """
        # Find the most significant non-zero chunk to avoid leading zeros in output
        first_digit_idx = -1
        for i in range(len(chunks) - 1, -1, -1):
            if chunks[i] != 0:
                first_digit_idx = i
                break
        
        if first_digit_idx == -1:
            return "0"

        # Convert the most significant chunk to string without padding
        s_num = str(chunks[first_digit_idx])
        
        # Convert the rest of the chunks, padding with zeros to maintain chunk size
        for i in range(first_digit_idx - 1, -1, -1):
            s_num += str(chunks[i]).zfill(chunk_size)
            
        return s_num

    # Step 3: The core multiplication logic, simulating the optimized C program
    def multiply_wuxing(p_chunks, q_chunks):
        """
        Multiplies two numbers represented as little-endian chunk lists.
        """
        o_chunks = [0] * O_CHUNKS

        for i in range(len(p_chunks)):
            carry = 0
            for j in range(len(q_chunks)):
                # In the Wuxing C program, 'product' would be a 'long' (10D)
                # to hold the result of multiplying two 'int's (5D).
                product = p_chunks[i] * q_chunks[j] + o_chunks[i+j] + carry
                o_chunks[i+j] = product % BASE
                carry = product // BASE
            
            # After the inner loop, propagate the final carry
            k = i + len(q_chunks)
            while carry > 0 and k < O_CHUNKS:
                sum_val = o_chunks[k] + carry
                o_chunks[k] = sum_val % BASE
                carry = sum_val // BASE
                k += 1
                
        return o_chunks

    # Step 4: Main execution flow
    # Define the input 100-digit numbers p and q
    p_str = "1234567890" * 10
    q_str = "9876543210" * 10

    # Convert input strings to the internal chunked representation
    p_arr = convert_str_to_chunks(p_str, P_DIGITS, CHUNK_DIGITS)
    q_arr = convert_str_to_chunks(q_str, Q_DIGITS, CHUNK_DIGITS)

    # Perform the multiplication
    o_arr = multiply_wuxing(p_arr, q_arr)

    # Convert the result back to a string for printing
    o_str = convert_chunks_to_str(o_arr, CHUNK_DIGITS)

    # Print the final equation as requested
    print(f"p = {p_str}")
    print(f"q = {q_str}")
    print(f"o = {o_str}")

solve_wuxing_multiplication()