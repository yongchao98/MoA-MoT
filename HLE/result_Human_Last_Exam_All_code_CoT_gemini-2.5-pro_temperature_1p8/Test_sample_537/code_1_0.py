def decode_word():
    """
    Decodes the given word based on the pattern:
    (number of unique consonants) - (number of unique vowels).
    """
    word = "anhaeketn"
    vowels = "aeiou"

    # Find the sets of unique characters, vowels, and consonants
    unique_chars = set(word)
    unique_vowels = {char for char in unique_chars if char in vowels}
    unique_consonants = unique_chars - unique_vowels

    # Get the counts
    num_unique_consonants = len(unique_consonants)
    num_unique_vowels = len(unique_vowels)

    # Calculate the final result
    result = num_unique_consonants - num_unique_vowels

    # Print the equation and the result
    print(f"{num_unique_consonants} - {num_unique_vowels} = {result}")

decode_word()