def find_chomsky_nouns():
    """
    This function identifies and prints the two nouns requested based on
    Noam Chomsky's examples in 'Syntactic Structures'.
    """
    print("In his 1957 work 'Syntactic Structures', Noam Chomsky lists several sentences to distinguish syntax from semantics.")
    print("The relevant list of sentences is as follows:")
    print("----------------------------------------------------------------")
    print("Grammatically correct sentences:")
    print("  1. Colorless green ideas sleep furiously.")
    print("  2. have you a book on modern music?")
    print("  3. the book seems interesting.")
    print("\nGrammatically incorrect sentences:")
    print("  4. *Furiously sleep ideas green colorless.")
    print("  5. *read you a book on modern music?")
    print("  6. *the book seems interestingly.")
    print("----------------------------------------------------------------\n")

    # The sentences in question
    last_correct_sentence = "the book seems interesting"
    last_incorrect_sentence = "the book seems interestingly"

    # The nouns from each sentence
    noun_from_correct = "book"
    noun_from_incorrect = "book"

    print(f"The last syntactically correct sentence is: '{last_correct_sentence}'")
    print(f"The noun in this sentence is: '{noun_from_correct}'")
    print("")
    print(f"The last syntactically incorrect sentence is: '{last_incorrect_sentence}'")
    print(f"The noun in this sentence is: '{noun_from_incorrect}'")

if __name__ == '__main__':
    find_chomsky_nouns()
