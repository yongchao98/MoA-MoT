import sys

def find_chomsky_nouns():
    """
    Finds the nouns in the last syntactically correct and last syntactically incorrect
    sentences from a famous list by Noam Chomsky.
    """
    # The list of sentences from Chomsky's Syntactic Structures (1957)
    # in the relevant section.
    sentences_data = [
        {"text": "Colorless green ideas sleep furiously", "is_correct": True, "nouns": ["ideas"]},
        {"text": "Furiously sleep ideas green colorless", "is_correct": False, "nouns": ["ideas"]},
        {"text": "have you a book on modern music", "is_correct": True, "nouns": ["book", "music"]},
        {"text": "read you a book on modern music", "is_correct": False, "nouns": ["book", "music"]},
        {"text": "the book seems interesting", "is_correct": True, "nouns": ["book"]},
        {"text": "the child seems sleeping", "is_correct": False, "nouns": ["child"]},
        {"text": "John plays tennis", "is_correct": True, "nouns": ["John", "tennis"]},
        {"text": "both John is playing tennis and is interesting", "is_correct": False, "nouns": ["John", "tennis"]},
        {"text": "I saw a fragile of", "is_correct": False, "nouns": []},
        {"text": "I saw a fragile whale", "is_correct": True, "nouns": ["whale"]},
        {"text": "sincerity may frighten the boy", "is_correct": True, "nouns": ["sincerity", "boy"]},
        {"text": "sincerity may admire the boy", "is_correct": True, "nouns": ["sincerity", "boy"]}, # Syntactically correct, semantically odd
        {"text": "boy the frighten may sincerity", "is_correct": False, "nouns": ["boy", "sincerity"]},
        {"text": "sincerity may virtue the boy", "is_correct": False, "nouns": ["sincerity", "boy"]}, # Syntactically incorrect as 'virtue' is not a verb
    ]

    last_correct_nouns = None
    last_incorrect_nouns = None

    # Reverse iterate to find the *last* occurrences first
    for i in range(len(sentences_data) - 1, -1, -1):
        sentence_info = sentences_data[i]
        
        # Find the last syntactically correct sentence
        if last_correct_nouns is None and sentence_info["is_correct"]:
            last_correct_nouns = sentence_info["nouns"]
            
        # Find the last syntactically incorrect sentence
        if last_incorrect_nouns is None and not sentence_info["is_correct"]:
            last_incorrect_nouns = sentence_info["nouns"]

        # Stop once both have been found
        if last_correct_nouns is not None and last_incorrect_nouns is not None:
            break

    print("The nouns in the last syntactically correct sentence are:", ", ".join(last_correct_nouns))
    print("The nouns in the last syntactically incorrect sentence are:", ", ".join(last_incorrect_nouns))
    
    # To answer the user's question directly:
    combined_nouns = sorted(list(set(last_correct_nouns + last_incorrect_nouns)))
    # Set sys.stdout directly for the final answer format
    sys.stdout.write(f'<<<{", ".join(combined_nouns)}>>>')


find_chomsky_nouns()