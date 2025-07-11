def find_chomsky_nouns():
    """
    Identifies the nouns in the last syntactically correct and incorrect
    sentences from the relevant section of Chomsky's "Syntactic Structures".
    """
    # In "Syntactic Structures" (1957), after introducing "Colorless green ideas sleep furiously,"
    # Chomsky provides a contrasting pair of sentences to further illustrate grammaticality.
    # These are the last two sentences in that specific introductory block of examples.

    last_syntactically_correct = "have you a book on modern music?"
    last_syntactically_incorrect = "read you a book on modern music?"

    # The nouns in these sentences are 'book' and 'music'. We can create a set
    # of these known nouns for identification.
    known_nouns = {"book", "music"}

    # We will process the words from both sentences to find the nouns used.
    # Using a set will automatically handle duplicates.
    words_from_sentences = last_syntactically_correct.split() + last_syntactically_incorrect.split()
    
    found_nouns = set()
    for word in words_from_sentences:
        # Clean the word of potential punctuation.
        cleaned_word = word.strip("?.,;")
        if cleaned_word in known_nouns:
            found_nouns.add(cleaned_word)
    
    # Sort the list for consistent output.
    final_nouns = sorted(list(found_nouns))

    print("The last syntactically correct sentence is: '{}'".format(last_syntactically_correct))
    print("The last syntactically incorrect sentence is: '{}'".format(last_syntactically_incorrect))
    print("\nThe two nouns used in these sentences are:")
    
    # Print each noun as a component of the final answer.
    if len(final_nouns) == 2:
        print(f"1: {final_nouns[0]}")
        print(f"2: {final_nouns[1]}")

if __name__ == '__main__':
    find_chomsky_nouns()