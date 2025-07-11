def vigenere_key(ciphertext, plaintext):
    """
    Calculates the Vigenère key used to encrypt a plaintext to a ciphertext.
    Key = (Ciphertext - Plaintext) mod 26.
    """
    key = []
    # Ensure inputs are lowercase for consistent alphabet mapping
    ciphertext = ciphertext.lower()
    plaintext = plaintext.lower()
    
    for i in range(len(ciphertext)):
        c_val = ord(ciphertext[i]) - ord('a')
        p_val = ord(plaintext[i]) - ord('a')
        # The formula for the key is (ciphertext - plaintext)
        k_val = (c_val - p_val + 26) % 26
        key.append(chr(k_val + ord('a')))
    return "".join(key)

def reverse_string(s):
    """Reverses a given string."""
    return s[::-1]

def find_p1():
    """
    Solves the recursive Vigenere cipher problem to find the original plaintext P_1.
    """
    # Given values from the problem statement
    p_1000 = "zuoeswzgnadou"
    e_1000 = "ikfcuwfgaoked"

    # From the problem definition, P_n = E_{n-1}. Therefore, P_1000 = E_999.
    e_999 = p_1000

    # We use a dictionary to store the E_n values as we compute them backward.
    E = {
        1000: e_1000,
        999: e_999
    }

    # The recursive relationship to find earlier encrypted strings is:
    # E_{n-2} = reverse(Vigenere_key(E_n, E_{n-1}))
    # We iterate from n = 1000 down to 3 to find all values down to E_1.
    for n in range(1000, 2, -1):
        e_n = E[n]
        e_n_minus_1 = E[n-1]
        
        # Calculate the key used at step n, which equals reverse(E_{n-2})
        key_n = vigenere_key(e_n, e_n_minus_1)
        
        # Reverse the key to get E_{n-2}
        e_n_minus_2 = reverse_string(key_n)
        
        # Store the computed value
        E[n-2] = e_n_minus_2
    
    # After the loop, we have E[2] and E[1]. We can now find P_1.
    # The relationship is: E_2 = Vigenere_encrypt(E_1, reverse(P_1))
    # This means the key for this operation is K_2 = reverse(P_1).
    # We find this key using: K_2 = Vigenere_key(E_2, E_1)
    e_2 = E[2]
    e_1 = E[1]
    
    key_2 = vigenere_key(e_2, e_1)
    
    # And since K_2 = reverse(P_1), we find P_1 by reversing K_2.
    p_1 = reverse_string(key_2)
    
    # Print the final result P_1 as requested.
    print(p_1)

# Run the function to find and print the answer.
find_p1()