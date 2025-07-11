def find_chomsky_nouns():
    """
    This function identifies and prints the nouns from the last syntactically correct
    and last syntactically incorrect sentences in the relevant section of
    Chomsky's "Syntactic Structures".
    """
    last_correct_sentence = "The child seems sleeping."
    last_incorrect_sentence = "The child seems sleep."

    # For these specific sentences, the nouns can be identified directly.
    # A more complex method like NLP part-of-speech tagging is not needed.
    noun_in_correct = "child"
    noun_in_incorrect = "child"

    print(f"The last syntactically correct sentence is: '{last_correct_sentence}'")
    print(f"The noun in this sentence is: '{noun_in_correct}'")
    print("-" * 20)
    print(f"The last syntactically incorrect sentence is: '{last_incorrect_sentence}'")
    print(f"The noun in this sentence is: '{noun_in_incorrect}'")

find_chomsky_nouns()