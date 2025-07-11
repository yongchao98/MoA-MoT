def solve_recursive_vigenere():
    """
    This script solves the recursive Vigenere cipher problem by working backward
    from the given P_1000 and E_1000 to find the original plaintext P_1.
    """

    # Helper function to reverse a string
    def reverse_string(s):
        return s[::-1]

    # Helper function to find the Vigenere key that encrypts a plaintext to a ciphertext.
    # The Vigenere encryption is C = (P + K) mod 26.
    # Therefore, the key is K = (C - P) mod 26.
    def find_key(ciphertext, plaintext):
        key = ""
        for i in range(len(ciphertext)):
            # Convert characters to numbers (a=0, b=1, ...)
            c_val = ord(ciphertext[i]) - ord('a')
            p_val = ord(plaintext[i]) - ord('a')
            
            # Calculate the key's character value
            k_val = (c_val - p_val + 26) % 26
            
            # Convert back to a character and append to the key
            key += chr(k_val + ord('a'))
        return key

    # --- Main Logic ---

    # Dictionary to store the calculated plaintexts P_n
    P = {}

    # Step 0: Use the given values
    P[1000] = "zuoeswzgnadou"
    e_1000 = "ikfcuwfgaoked"

    # Step 1: Calculate P[999]
    # The relation is E_1000 = Encrypt(P_1000, reverse(P_999)).
    # So, we find the key (which is reverse(P_999)) that turns P_1000 into E_1000.
    key_1000 = find_key(e_1000, P[1000])
    # Then we reverse the key to get P_999.
    P[999] = reverse_string(key_1000)

    # Step 2: Iterate backwards from n=1000 down to n=3
    # In each step, we calculate P[n-2] from P[n] and P[n-1].
    # The relation is P[n-1] = Decrypt(P[n], reverse(P[n-2])).
    # So, reverse(P[n-2]) is the key that turns P[n] into P[n-1].
    for n in range(1000, 2, -1):
        key_n_minus_1 = find_key(P[n], P[n-1])
        P[n-2] = reverse_string(key_n_minus_1)

    # Step 3: Present the final calculation for P_1
    p_3 = P[3]
    p_2 = P[2]
    p_1 = P[1]

    print("To find P_1, we work backward recursively.")
    print("The final step in this process is calculating P_1 from P_3 and P_2.")
    print("-" * 50)
    print(f"Calculated P_3 = \"{p_3}\"")
    print(f"Calculated P_2 = \"{p_2}\"")
    print("-" * 50)
    print("The relationship is P_3 = Encrypt(P_2, reverse(P_1)).")
    print("To find reverse(P_1), we find the key (K_2) that encrypts P_2 to P_3.")
    print("Using the formula: Key = (Ciphertext - Plaintext) mod 26.\n")
    print("Let's calculate K_2 = reverse(P_1) character by character:")
    
    key_2 = find_key(p_3, p_2)
    
    for i in range(len(p_3)):
        p3_char = p_3[i]
        p2_char = p_2[i]
        k2_char = key_2[i]
        p3_val = ord(p3_char) - ord('a')
        p2_val = ord(p2_char) - ord('a')
        k2_val = (p3_val - p2_val + 26) % 26
        print(f"K_2[{i:02d}]: ('{p3_char}' - '{p2_char}') mod 26 = ({p3_val:2d} - {p2_val:2d}) mod 26 = {k2_val:2d} => '{k2_char}'")

    print(f"\nThis gives the key K_2 = reverse(P_1) = \"{key_2}\"")
    print(f"Finally, we reverse K_2 to find P_1: P_1 = reverse(\"{key_2}\")")
    print(f"\nThe final answer is P_1 = \"{p_1}\"")
    
    # Final answer in the required format
    print(f"\n<<<" + p_1 + ">>>")

solve_recursive_vigenere()