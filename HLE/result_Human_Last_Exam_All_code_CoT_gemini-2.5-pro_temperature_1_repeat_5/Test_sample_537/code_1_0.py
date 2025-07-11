def decode_word(word):
    """
    Decodes the word by counting its consonants and prints the equation.
    """
    vowels = "aeiou"
    
    # Lists to store consonants and the numbers for the equation
    consonants_found = []
    equation_numbers = []

    # Iterate over each character in the word to find consonants
    for char in word:
        if char.isalpha() and char.lower() not in vowels:
            consonants_found.append(char)
            equation_numbers.append("1")

    # The decoded number is the total count of consonants
    decoded_number = len(consonants_found)

    # Building the equation string as requested
    equation_str = " + ".join(equation_numbers)
    
    print(f"Decoding the word: '{word}'")
    print(f"The consonants are: {', '.join(consonants_found)}")
    print("The final equation is the sum of 1 for each consonant:")
    print(f"{equation_str} = {decoded_number}")

# The word to be decoded
last_word = "anhaeketn"

# Run the decoding process
decode_word(last_word)