import textwrap

def wuxing_rsa_multiply():
    """
    This program simulates an optimized multiplication of two 100-digit numbers
    on the conceptual Wuxing computer architecture.
    """
    # --- Wuxing Architecture & Problem Constants ---
    INT_SIZE_D = 5
    INT_BASE = 10**INT_SIZE_D  # 100000, our calculation base
    P_DIGITS = 100
    Q_DIGITS = 100
    P_CHUNKS = P_DIGITS // INT_SIZE_D
    Q_CHUNKS = Q_DIGITS // INT_SIZE_D
    O_CHUNKS = (P_DIGITS + Q_DIGITS) // INT_SIZE_D

    # --- Input Numbers ---
    # Two 100-digit numbers, p and q.
    p_str = "1212121212121212121212121212121212121212121212121212121212121212121212121212121212121212121212121212"
    q_str = "8787878787878787878787878787878787878787878787878787878787878787878787878787878787878787878787878787"

    # --- Step 1: Parse input strings into int arrays (chunks) ---
    # This simulates creating int p_arr[20] and int q_arr[20] in C.
    p_arr = [0] * P_CHUNKS
    q_arr = [0] * Q_CHUNKS

    for i in range(P_CHUNKS):
        start = max(0, P_DIGITS - (i + 1) * INT_SIZE_D)
        end = P_DIGITS - i * INT_SIZE_D
        p_arr[i] = int(p_str[start:end])

    for i in range(Q_CHUNKS):
        start = max(0, Q_DIGITS - (i + 1) * INT_SIZE_D)
        end = Q_DIGITS - i * INT_SIZE_D
        q_arr[i] = int(q_str[start:end])

    # --- Step 2: Perform Optimized Multiplication ---
    # This simulates the grade-school algorithm:
    # int o_arr[40] = {0};
    # for (j=0; j<20; j++) {
    #   long carry = 0;
    #   for (i=0; i<20; i++) {
    #     long product = (long)p_arr[i] * q_arr[j] + o_arr[i+j] + carry;
    #     o_arr[i+j] = product % INT_BASE;
    #     carry = product / INT_BASE;
    #   }
    #   o_arr[j+20] += carry;
    # }
    o_arr = [0] * O_CHUNKS

    for j in range(Q_CHUNKS):
        carry = 0  # Simulates a 'long' register for the carry
        for i in range(P_CHUNKS):
            # product simulates a 'long' register calculation. Python handles large integers.
            product = p_arr[i] * q_arr[j] + o_arr[i + j] + carry
            o_arr[i + j] = product % INT_BASE
            carry = product // INT_BASE
        
        # After iterating through p_arr, propagate the final carry
        if j + P_CHUNKS < O_CHUNKS:
            o_arr[j + P_CHUNKS] += carry


    # --- Step 3: Format the result for output ---
    output_str = ""
    # Find the most significant non-zero chunk
    first_chunk_idx = O_CHUNKS - 1
    while first_chunk_idx >= 0 and o_arr[first_chunk_idx] == 0:
        first_chunk_idx -= 1

    if first_chunk_idx == -1:
        output_str = "0"
    else:
        # Append the first chunk without padding
        output_str += str(o_arr[first_chunk_idx])
        # Append subsequent chunks with 5-digit zero-padding
        for i in range(first_chunk_idx - 1, -1, -1):
            output_str += f"{o_arr[i]:05d}"
            
    # --- Print Final Equation ---
    print("p = " + p_str)
    print("q = " + q_str)
    # Using textwrap to make the large output readable
    print("o = p * q = " + "\n".join(textwrap.wrap(output_str, 80)))


# Execute the simulation
wuxing_rsa_multiply()
