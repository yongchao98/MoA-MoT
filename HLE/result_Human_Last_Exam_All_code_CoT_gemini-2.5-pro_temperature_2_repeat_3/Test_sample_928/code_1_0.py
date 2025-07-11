import nltk

# The NLTK library is used for natural language processing.
# If you don't have it, you can install it (`pip install nltk`)
# and download the necessary data packages by running the following
# two lines of code in a Python interpreter once:
# nltk.download('punkt')
# nltk.download('averaged_perceptron_tagger')

def find_nouns_in_chomsky_sentences():
    """
    Identifies the nouns in the last correct and incorrect sentences
    Chomsky uses after "Colorless green ideas sleep furiously" in
    Syntactic Structures (1957).
    """
    # In Syntactic Structures, after the "Colorless green ideas" example,
    # Chomsky provides these two sentences to further his point.
    last_correct_sentence = "Have you a book on modern music"
    last_incorrect_sentence = "Read you a book on modern music"

    # Combine text for analysis
    text = last_correct_sentence + ". " + last_incorrect_sentence + "."

    # Tokenize the combined text into words
    tokens = nltk.word_tokenize(text.lower())

    # Use Part-of-Speech (POS) tagging to identify the role of each word
    tagged_words = nltk.pos_tag(tokens)

    # A set is used to store the unique nouns found.
    # We are looking for singular nouns, which are tagged as 'NN'.
    found_nouns = set()
    for word, tag in tagged_words:
        if tag == 'NN':
            found_nouns.add(word)

    # Print the results clearly. The list is sorted for consistent output.
    sorted_nouns = sorted(list(found_nouns))
    print("The two nouns from the last syntactically correct sentence ('Have you a book on modern music?')")
    print("and the last syntactically incorrect sentence ('Read you a book on modern music?') are:")
    print(f"{sorted_nouns[0]}")
    print(f"{sorted_nouns[1]}")

find_nouns_in_chomsky_sentences()