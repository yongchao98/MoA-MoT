def solve_cipher_puzzle():
    """
    Solves the recursive Vigenere cipher problem by working backward
    from the 1000th step to find the original plaintext P_1.
    """

    def vigenere_find_key(ciphertext, plaintext):
        """
        Calculates the Vigenere key given a ciphertext and plaintext.
        The formula is K_i = (E_i - P_i) mod 26.
        """
        key = []
        for i in range(len(ciphertext)):
            e_val = ord(ciphertext[i]) - ord('a')
            p_val = ord(plaintext[i]) - ord('a')
            k_val = (e_val - p_val + 26) % 26
            key.append(chr(k_val + ord('a')))
        return "".join(key)

    # Initial known values at step n=1000
    p_1000 = "zuoeswzgnadou"
    e_1000 = "ikfcuwfgaoked"

    # Initialize variables for the backward iteration.
    p_current = p_1000
    e_current = e_1000

    # Iterate backwards from n = 1000 down to 2.
    # In each step 'n', we calculate the values for step 'n-1'.
    for n in range(1000, 1, -1):
        # 1. Find K_n using the current P_n and E_n.
        k_n = vigenere_find_key(e_current, p_current)
        
        # 2. From the problem, K_n = reverse(P_{n-1}), so P_{n-1} = reverse(K_n).
        p_previous = k_n[::-1]
        
        # 3. From the problem, E_{n-1} = P_n.
        e_previous = p_current
        
        # 4. Update state for the next iteration (which will be for step n-1).
        p_current = p_previous
        e_current = e_previous

    # After the loop, p_current holds the value for P_1.
    p_1 = p_current

    # Per the instructions, convert the final answer string into numbers
    # (a=0, b=1, ...) and print them.
    final_numbers = [ord(char) - ord('a') for char in p_1]
    
    # The final equation is P_1 = <result>. The numbers in this equation
    # are the numeric values of the characters in the result.
    print(*final_numbers)

solve_cipher_puzzle()