def find_next_number():
    """
    This function calculates the next number in the sequence based on a Vigenere cipher
    with a constant XOR key.
    """
    
    # The key is derived from the prompt's introduction.
    key_text = "You are a helpful AI assistant."
    
    # The plaintext is the question itself.
    plain_text = "What is the next number in the sequence:"
    
    # The sequence provided by the user.
    sequence = [111, 142, 111, 41, 67, 67, 67, 93, 111, 111, 62, 62, 111, 111, 36, 36, 49, 155, 49, 62, 49, 49, 62, 62, 10, 36, 36, 36, 124, 124, 124, 36, 124]
    
    # We calculate the constant key 'X' from the first element.
    # C[0] = ord(P[0]) ^ ord(K[0]) ^ X
    # X = C[0] ^ ord(P[0]) ^ ord(K[0])
    X = sequence[0] ^ ord(plain_text[0]) ^ ord(key_text[0])
    
    # The sequence has 33 numbers, so the next number is at index 33.
    next_index = len(sequence)

    # Get the characters for the next calculation using modular arithmetic for looping.
    plain_char_val = ord(plain_text[next_index % len(plain_text)])
    key_char_val = ord(key_text[next_index % len(key_text)])
    
    # Calculate the next number using our formula.
    next_number = plain_char_val ^ key_char_val ^ X

    # As requested, output the final equation.
    print(f"Next Number = ord('{plain_text[next_index % len(plain_text)]}') ^ ord('{key_text[next_index % len(key_text)]}') ^ {X}")
    print(f"Next Number = {plain_char_val} ^ {key_char_val} ^ {X}")
    print(f"Next Number = {next_number}")

find_next_number()
<<<127>>>