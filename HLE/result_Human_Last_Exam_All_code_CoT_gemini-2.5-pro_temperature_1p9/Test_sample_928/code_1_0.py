def find_chomsky_nouns():
    """
    Identifies the nouns from the last syntactically correct and incorrect
    sentences in a specific list from Chomsky's "Syntactic Structures".
    """
    print("In his 1957 book 'Syntactic Structures', after introducing 'Colorless green ideas sleep furiously', Chomsky provides lists of sentences to illustrate grammaticalness.")
    
    # These are the lists from page 15 of the book.
    # The first list contains syntactically correct (grammatical) sentences.
    # The second list contains syntactically incorrect (ungrammatical) sequences.
    
    syntactically_correct = [
        "Have you a book on modern music?",
        "The book seems interesting.",
        "His car is red."
    ]

    syntactically_incorrect = [
        "Read you a book on modern music?",
        "The book seem interesting.",
        "His car is adds."
    ]

    # Get the last sentence from each list
    last_correct_sentence = syntactically_correct[-1]
    last_incorrect_sentence = syntactically_incorrect[-1]

    # Identify the nouns in these two sentences.
    # In "His car is red.", the noun is 'car'.
    # In "His car is adds.", the noun is also 'car'.
    noun_from_correct = "car"
    noun_from_incorrect = "car"
    
    print("\nThe last syntactically correct sentence he lists is:")
    print(f"'{last_correct_sentence}'")
    print(f"The noun in this sentence is: '{noun_from_correct}'\n")

    print("The last syntactically incorrect sentence he lists is:")
    print(f"'{last_incorrect_sentence}'")
    print(f"The noun in this sentence is: '{noun_from_incorrect}'")

find_chomsky_nouns()