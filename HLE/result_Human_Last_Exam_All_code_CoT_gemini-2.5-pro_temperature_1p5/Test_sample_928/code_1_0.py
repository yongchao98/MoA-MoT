def find_chomsky_nouns():
    """
    This script identifies the two nouns used in the last syntactically correct
    and last syntactically incorrect sentences from a famous list by Noam Chomsky.
    """

    # In "Syntactic Structures" (1957), Chomsky provides two lists of sentences.
    # The following are the last sentences from each list.
    last_syntactically_correct = "sincerity may frighten the boy"
    last_syntactically_incorrect = "the boy may frighten sincerity"

    # The nouns in the first sentence are 'sincerity' and 'boy'.
    noun1_sentence1 = "sincerity"
    noun2_sentence1 = "boy"

    # The nouns in the second sentence are 'boy' and 'sincerity'.
    noun1_sentence2 = "boy"
    noun2_sentence2 = "sincerity"

    print(f"The last syntactically correct sentence is: '{last_syntactically_correct}'")
    print(f"The nouns in this sentence are: '{noun1_sentence1}' and '{noun2_sentence1}'")
    print("-" * 20)
    print(f"The last syntactically incorrect sentence is: '{last_syntactically_incorrect}'")
    print(f"The nouns in this sentence are: '{noun1_sentence2}' and '{noun2_sentence2}'")
    print("-" * 20)
    print(f"\nThe two nouns used in these two sentences are '{noun1_sentence1}' and '{noun2_sentence1}'.")

find_chomsky_nouns()