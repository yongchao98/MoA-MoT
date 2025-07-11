def find_chomsky_nouns():
    """
    Identifies and prints the two nouns from the last syntactically correct
    and last syntactically incorrect sentences Chomsky gives in the section
    of "Syntactic Structures" where "Colorless green ideas sleep furiously" is introduced.
    """
    # In "Syntactic Structures", Chomsky lists several grammatical but semantically odd sentences.
    # The last one often cited in this context is:
    last_correct_sentence = "sincerity may admire the boy"

    # He contrasts these with ungrammatical sentences. A corresponding ungrammatical version is:
    last_incorrect_sentence = "sincerity may the boy admire"

    # The nouns in both of these sentences are the same.
    noun1 = "boy"
    noun2 = "sincerity"

    print("In his work 'Syntactic Structures', Chomsky provides lists of sentences to distinguish syntax from semantics.")
    print("-" * 80)
    print(f"The last syntactically CORRECT sentence in the list is often cited as:")
    print(f'"{last_correct_sentence}"')
    print("\n")
    print(f"The corresponding syntactically INCORRECT sentence is:")
    print(f'"{last_incorrect_sentence}"')
    print("-" * 80)
    print(f"The two nouns used in these sentences are: '{noun1}' and '{noun2}'.")

find_chomsky_nouns()