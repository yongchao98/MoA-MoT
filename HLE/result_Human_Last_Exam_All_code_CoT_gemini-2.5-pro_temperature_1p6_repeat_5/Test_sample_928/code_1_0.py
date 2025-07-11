def find_chomsky_nouns():
    """
    Identifies nouns in specific sentences from Chomsky's work.

    In "Syntactic Structures", Chomsky provides several examples to distinguish
    syntax from semantics. A representative list, based on the examples he gives, is:
    1. Colorless green ideas sleep furiously. (Correct)
    2. Sincerity may frighten the boy. (Correct)
    3. *Furiously sleep ideas green colorless. (Incorrect)
    4. *Sincerity may virtue the boy. (Incorrect)

    This code identifies the nouns in the last correct and last incorrect
    sentences from such a list.
    """

    # The last syntactically correct sentence in the relevant set of examples.
    last_correct_sentence = "Sincerity may frighten the boy."
    # The last syntactically incorrect sentence.
    last_incorrect_sentence = "*Sincerity may virtue the boy."

    # Nouns found in "Sincerity may frighten the boy".
    nouns_in_correct_sentence = ["Sincerity", "boy"]

    # Nouns found in "*Sincerity may virtue the boy".
    # Here, "virtue" is used ungrammatically as a verb.
    nouns_in_incorrect_sentence = ["Sincerity", "boy"]

    # The problem asks for the nouns used in both sentences.
    # We can use a set to find the unique nouns across both lists.
    combined_nouns = set(nouns_in_correct_sentence + nouns_in_incorrect_sentence)

    print(f"The last syntactically correct sentence is: '{last_correct_sentence}'")
    print(f"The last syntactically incorrect sentence is: '{last_incorrect_sentence}'")
    print("\nThe two nouns used in these sentences are:")

    # Print each of the final nouns.
    for noun in sorted(list(combined_nouns)):
        print(noun)

find_chomsky_nouns()