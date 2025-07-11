def vigenere_get_key(plaintext, ciphertext):
    """
    Calculates the Vigenere key given a plaintext and a corresponding ciphertext.
    The formula for each character is: K_i = (C_i - P_i) mod 26.
    """
    key = ""
    for p_char, c_char in zip(plaintext, ciphertext):
        # Convert characters to numbers (a=0, z=25)
        p_num = ord(p_char) - ord('a')
        c_num = ord(c_char) - ord('a')
        
        # Calculate the key character's number
        k_num = (c_num - p_num) % 26
        
        # Convert number back to character and append to the key
        key += chr(ord('a') + k_num)
    return key

def solve_crypto_puzzle():
    """
    Works backwards from P_1000 and E_1000 to find the original P_1.
    """
    # Given values at the end of the chain
    p_1000 = "zuoeswzgnadou"
    e_1000 = "ikfcuwfgaoked"

    # --- Step 1: Find P_999 ---
    # The relation is: E_1000 = Encrypt(P_1000, reverse(P_999))
    # We find the key used, which is reverse(P_999).
    key_for_999 = vigenere_get_key(p_1000, e_1000)
    # Then we reverse the key to get P_999.
    p_999 = key_for_999[::-1]
    
    # Initialize variables for the loop. We will work backwards from n=1000.
    # At the start of each loop iteration 'n', p_n will hold P_n and p_n_minus_1 will hold P_{n-1}.
    p_n = p_1000
    p_n_minus_1 = p_999

    # --- Step 2: Iterate backwards from n=1000 down to n=3 ---
    # The general backward recurrence relation is: P_{n-2} = reverse(find_key(P_{n-1}, P_n))
    for n in range(1000, 2, -1):
        # In this loop for n, we calculate P_{n-2}
        key = vigenere_get_key(p_n_minus_1, p_n)
        p_n_minus_2 = key[::-1]
        
        # Update variables for the next iteration (i.e., for n-1)
        p_n = p_n_minus_1
        p_n_minus_1 = p_n_minus_2
        
    # After the loop finishes (last run was for n=3), p_n_minus_1 will hold the calculated P_1.
    p_1 = p_n_minus_1
    
    # Print the final result. No equation is involved, so we print the calculated string.
    print(f"P_1 = {p_1}")

solve_crypto_puzzle()
<<<iamthevigenere>>>