def find_key(ciphertext, plaintext):
    """
    Calculates the Vigenere key given a ciphertext and plaintext.
    Assumes all characters are lowercase English letters.
    """
    key = []
    # Ensure the key is the same length as the plaintext/ciphertext
    for i in range(len(plaintext)):
        # Convert characters to numbers (a=0, b=1, ...)
        cipher_val = ord(ciphertext[i]) - ord('a')
        plain_val = ord(plaintext[i]) - ord('a')
        
        # Calculate key character value using (C - P) mod 26
        # Add 26 to handle negative results before the modulo operation
        key_val = (cipher_val - plain_val + 26) % 26
        
        # Convert number back to character and append to the key
        key.append(chr(key_val + ord('a')))
        
    return "".join(key)

def solve_cipher():
    """
    Works backward from P_1000 and E_1000 to find P_1.
    """
    # Initial values for n=1000
    p_n = "zuoeswzgnadou" # This is P_1000
    e_n = "ikfcuwfgaoked" # This is E_1000

    # We need to find P_1. We iterate backward from n=1000 down to n=2.
    # In each iteration, we find P_{n-1}.
    for n in range(1000, 1, -1):
        # At any step n (for n>1), the relationship is:
        # E_n = VigenereEncrypt(P_n, K_n) where K_n = reverse(P_{n-1})
        
        # We can find K_n from E_n and P_n.
        key_n = find_key(e_n, p_n)
        
        # We know K_n is the reversed version of P_{n-1}.
        # So, to get P_{n-1}, we reverse K_n.
        p_n_minus_1 = key_n[::-1]
        
        # Now we set up the values for the next iteration (n-1).
        # The new plaintext is P_{n-1}.
        # The new ciphertext E_{n-1} is equal to the old plaintext P_n.
        p_n = p_n_minus_1
        e_n = p_n
        # The loop continues until p_n holds the value of P_1.
    
    # After the loop, p_n contains P_1
    print("The original plaintext P_1 is:")
    print(p_n)

# Run the solver
solve_cipher()