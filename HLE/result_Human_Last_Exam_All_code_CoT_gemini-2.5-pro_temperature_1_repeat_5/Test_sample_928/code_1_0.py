def find_chomsky_nouns():
    """
    This script identifies the two nouns from the last syntactically correct sentence
    and the last syntactically incorrect sentence in the list Chomsky provides
    in 'Syntactic Structures' to illustrate the syntax vs. semantics distinction.
    """

    # The list of example sentences from Chomsky's 'Syntactic Structures'
    # used to demonstrate different levels of grammaticalness.
    correct_sentences = [
        "sincerity may frighten the boy",
        "the boy may admire sincerity"
    ]

    incorrect_sentences = [
        "boy the may frighten sincerity",
        "sincerity may admire boy the"
    ]

    # The last syntactically correct sentence given in his example list
    last_correct_sentence = correct_sentences[-1]

    # The last syntactically incorrect sentence given in his example list
    last_incorrect_sentence = incorrect_sentences[-1]

    # The nouns are the same in both of these sentences
    noun1 = "sincerity"
    noun2 = "boy"

    print("In the section of 'Syntactic Structures' discussing syntax vs. semantics, Chomsky provides several examples.")
    print(f"The last syntactically correct sentence he lists is: '{last_correct_sentence}'")
    print(f"The last syntactically incorrect sentence he lists is: '{last_incorrect_sentence}'")
    print("\nThe two nouns used in these sentences are:")
    print(f"1. {noun1}")
    print(f"2. {noun2}")

find_chomsky_nouns()