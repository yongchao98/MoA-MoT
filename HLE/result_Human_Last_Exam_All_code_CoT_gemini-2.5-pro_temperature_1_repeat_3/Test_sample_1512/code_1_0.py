import math

def multiply_large_numbers(p_str, q_str):
    """
    Multiplies two large numbers represented as strings, using an algorithm
    optimized for the hypothetical Wuxing architecture. This function simulates
    multiplying the numbers in 5-digit chunks.
    """
    # On Wuxing, we'd use 5D 'int' types.
    # The base for our chunks is 10^5 = 100000.
    BASE = 100000
    CHUNK_SIZE = 5

    # 1. Parse input strings into lists of 5-digit integer chunks.
    # The numbers are stored with the least significant chunk at index 0 for easier processing.
    p_chunks = []
    for i in range(len(p_str), 0, -CHUNK_SIZE):
        start = max(0, i - CHUNK_SIZE)
        p_chunks.append(int(p_str[start:i]))

    q_chunks = []
    for i in range(len(q_str), 0, -CHUNK_SIZE):
        start = max(0, i - CHUNK_SIZE)
        q_chunks.append(int(q_str[start:i]))
    
    # 2. Initialize the result list.
    # A 100-digit * 100-digit number can have up to 200 digits.
    # Number of chunks for p: ceil(100/5) = 20
    # Number of chunks for q: ceil(100/5) = 20
    # Number of chunks for the result: 20 + 20 = 40
    len_p = len(p_chunks)
    len_q = len(q_chunks)
    o_chunks = [0] * (len_p + len_q)

    # 3. Perform multiplication using the chunk-based "schoolbook" method.
    for i in range(len_q):
        carry = 0
        for j in range(len_p):
            # In Wuxing, this product would be calculated using a 10D 'long'.
            # The calculation is: (p_chunk * q_chunk) + existing_result_chunk + carry
            product = p_chunks[j] * q_chunks[i] + o_chunks[i+j] + carry
            o_chunks[i+j] = product % BASE
            carry = product // BASE
        
        # After the inner loop, place the final carry in the next position.
        if i + len_p < len(o_chunks):
            o_chunks[i + len_p] += carry

    # 4. Format the result list back into a final string.
    # Remove any leading zero chunks from the result (e.g., if inputs were small).
    while len(o_chunks) > 1 and o_chunks[-1] == 0:
        o_chunks.pop()

    # The most significant chunk is converted to a string without padding.
    result_str = str(o_chunks[-1])
    
    # All other chunks are formatted to 5 digits, padded with leading zeros if necessary.
    for i in range(len(o_chunks) - 2, -1, -1):
        result_str += f"{o_chunks[i]:0{CHUNK_SIZE}d}"
        
    return result_str

# Define the 100-digit input numbers p and q.
p_input = "1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890"
q_input = "9876543210987654321098765432109876543210987654321098765432109876543210987654321098765432109876543210"

# Calculate the product using our Wuxing-optimized method.
o_output = multiply_large_numbers(p_input, q_input)

# Print the final equation with all numbers.
print(f"{p_input} * {q_input} = {o_output}")