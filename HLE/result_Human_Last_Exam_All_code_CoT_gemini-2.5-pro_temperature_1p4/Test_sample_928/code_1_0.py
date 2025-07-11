def find_chomsky_nouns():
    """
    Identifies and prints the two nouns from Chomsky's specified sentences.
    
    In his 1957 book "Syntactic Structures", Noam Chomsky provides a list of 
    sentences to demonstrate the independence of syntax from semantics.
    The sentences are:
    1. "Colorless green ideas sleep furiously." (Correct)
    2. "Furiously sleep ideas green colorless." (Incorrect)
    3. "have you a book on modern music?" (Correct)
    4. "the book seems interesting." (Correct)
    5. "I saw a fragile of." (Incorrect)
    6. "I saw a fragile whale." (Correct)
    """

    # The last syntactically correct sentence listed is #6.
    last_correct_sentence = "I saw a fragile whale."
    noun_from_correct = "whale"

    # The incorrect sentences are #2 and #5. The last one (#5) has no noun.
    # Therefore, we take the noun from the other incorrect sentence, #2.
    incorrect_sentence_with_noun = "Furiously sleep ideas green colorless."
    noun_from_incorrect = "ideas"

    print("Identifying the nouns based on Noam Chomsky's sentences in 'Syntactic Structures':")
    print("-" * 70)
    
    print(f"The last syntactically correct sentence is: '{last_correct_sentence}'")
    print(f"The noun from this sentence is: '{noun_from_correct}'")
    print("")
    
    print("The last syntactically incorrect sentence is 'I saw a fragile of', which contains no nouns.")
    print(f"The other incorrect sentence listed is: '{incorrect_sentence_with_noun}'")
    print(f"The noun from this sentence is: '{noun_from_incorrect}'")
    print("-" * 70)

    print(f"\nThe two requested nouns are '{noun_from_incorrect}' and '{noun_from_correct}'.")

find_chomsky_nouns()