def decode_words():
    """
    Finds the pattern in a list of words and decodes them.
    The pattern is that the decoded word is the English name
    for the number of letters in the encoded word.
    """
    words = [
        "utpshtheas",
        "fkuhu",
        "numhkatnatu",
        "anhaeketn"
    ]

    # A mapping from number to its English word representation
    num_to_word_map = {
        1: "ONE", 2: "TWO", 3: "THREE", 4: "FOUR", 5: "FIVE",
        6: "SIX", 7: "SEVEN", 8: "EIGHT", 9: "NINE", 10: "TEN",
        11: "ELEVEN", 12: "TWELVE"
    }

    print("Decoding based on the pattern found:")
    for word in words:
        length = len(word)
        decoded_word = num_to_word_map.get(length, "UNKNOWN")
        print(f"'{word}' has a length of {length}, so it decodes to '{decoded_word}'")
    
    # Specifically for the last word, showing the 'equation' as requested.
    last_word = words[-1]
    last_word_length = len(last_word)
    final_decoded_word = num_to_word_map.get(last_word_length, "UNKNOWN")
    
    print("\n--- Final Answer ---")
    print("The final equation is: length('{}') = {}".format(last_word, last_word_length))
    print("The decoded word is: {}".format(final_decoded_word))

if __name__ == '__main__':
    decode_words()