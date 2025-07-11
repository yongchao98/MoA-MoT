def find_chomsky_nouns():
    """
    This script identifies the two nouns from the last pair of sentences
    Chomsky used in the section of "Syntactic Structures" where he
    introduced "Colorless green ideas sleep furiously."
    """

    last_correct_sentence = "Have you a book on modern music?"
    last_incorrect_sentence = "Read you a book on modern music?"

    # The nouns are identical in both sentences.
    # In a real-world scenario, a natural language processing library (like NLTK or SpaCy)
    # would be used for part-of-speech tagging. For this specific problem,
    # we can identify them directly.
    nouns_in_sentences = {"book", "music"}

    noun1 = "book"
    noun2 = "music"

    print(f"The last syntactically correct sentence is: '{last_correct_sentence}'")
    print(f"The last syntactically incorrect sentence is: '{last_incorrect_sentence}'")
    print("\nThe two nouns used in these sentences are:")
    print(noun1)
    print(noun2)


find_chomsky_nouns()