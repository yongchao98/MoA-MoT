def solve():
    """
    Decodes the given ciphertext using an Autokey cipher.
    The formula is P[i] = (C[i] - P[i-1]) mod 26.
    The key P[0] is derived from the first letter of the ciphertext.
    """
    ciphertext = "anhaeketn"
    
    # Step 1: Determine the key P[0]
    # The key is the first letter of the ciphertext.
    key_char = ciphertext[0]
    key_val = ord(key_char) - ord('a')
    
    # Initialize the list of plaintext character values
    p_vals = []
    
    # The previous plaintext character value, starts with the key
    prev_p_val = key_val
    
    print(f"Decoding '{ciphertext}' with key '{key_char}' ({key_val})")
    print("-" * 30)

    # Step 2: Decode each character
    for i, c_char in enumerate(ciphertext):
        c_val = ord(c_char) - ord('a')
        
        # Apply the formula: P[i] = C[i] - P[i-1]
        p_val = (c_val - prev_p_val + 26) % 26
        p_char = chr(p_val + ord('a'))
        
        # Print the equation for the current character
        print(f"'{p_char}' = ('{c_char}' - '{chr(prev_p_val + ord('a'))}') mod 26")
        print(f"  {p_val} = ({c_val} - {prev_p_val}) mod 26")
        
        # Store the result and update the previous plaintext value for the next iteration
        p_vals.append(p_val)
        prev_p_val = p_val

    # Step 3: Construct the final decoded word
    decoded_word = "".join([chr(val + ord('a')) for val in p_vals])
    
    print("-" * 30)
    print(f"The decoded word is: {decoded_word}")

solve()
<<<vigenere>>>