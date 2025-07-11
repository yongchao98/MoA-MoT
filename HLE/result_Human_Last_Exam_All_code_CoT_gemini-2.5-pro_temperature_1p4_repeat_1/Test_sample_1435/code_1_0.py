def vigenere_subtract(p_str, e_str):
    """
    Calculates the key K from plaintext P and ciphertext E.
    The Vigenère relation is E = P + K (mod 26).
    Therefore, K = E - P (mod 26).
    """
    result = []
    for p_char, e_char in zip(p_str, e_str):
        p_val = ord(p_char) - ord('a')
        e_val = ord(e_char) - ord('a')
        k_val = (e_val - p_val + 26) % 26
        result.append(chr(k_val + ord('a')))
    return "".join(result)

def reverse_string(s):
    """Reverses a string."""
    return s[::-1]

def solve_for_p1():
    """
    Solves the recursive Vigenère cipher problem to find P_1.
    """
    # Given values at step 1000
    p1000 = "zuoeswzgnadou"
    e1000 = "ikfcuwfgaoked"

    # To start, we calculate P_999.
    # The relation is E_1000 = Encrypt(P_1000, reverse(P_999)).
    # So, reverse(P_999) = Subtract(P_1000, E_1000).
    key_for_1000 = vigenere_subtract(p1000, e1000)
    p999 = reverse_string(key_for_1000)

    # We now have a recurrence relation P_{i-2} = reverse(Subtract(P_{i-1}, P_i)).
    # We will use two variables to track the last two plaintexts in the sequence.
    # Let p_i be P_i and p_i_minus_1 be P_{i-1}.
    p_i = p1000
    p_i_minus_1 = p999

    # We need to find P_1. We have P_1000 and P_999.
    # We loop from i=1000 down to 3 to find P_{i-2}. This is 998 iterations.
    for _ in range(998):
        # Calculate P_{i-2} using P_{i-1} and P_i
        # reverse(P_{i-2}) = Subtract(P_{i-1}, P_i)
        key_i_minus_1 = vigenere_subtract(p_i_minus_1, p_i)
        p_i_minus_2 = reverse_string(key_i_minus_1)
        
        # Update variables for the next iteration.
        # The new P_i becomes the current P_{i-1}.
        # The new P_{i-1} becomes the calculated P_{i-2}.
        p_i = p_i_minus_1
        p_i_minus_1 = p_i_minus_2
    
    # After the loop, p_i_minus_1 will hold the value of P_1.
    p1 = p_i_minus_1
    print(p1)

solve_for_p1()