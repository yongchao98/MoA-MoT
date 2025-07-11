import collections

def find_similar_words():
    """
    Identifies and prints a comparison of words from two Asian languages
    with no direct historical contact.
    """
    # The two languages chosen are Korean and Tamil.
    # - They are from different, unrelated language families (Koreanic and Dravidian).
    # - They are geographically distant and have no significant shared history of contact.
    language1 = "Korean"
    language2 = "Tamil"

    # Using an ordered dictionary to store word pairs for a clean printout.
    word_pairs = collections.OrderedDict()
    word_pairs["Mom"] = {"language": language1, "word": "eomma (엄마)", "similarity_with": language2, "similar_word": "ammā (அம்மா)"}
    word_pairs["Dad"] = {"language": language1, "word": "appa (아빠)", "similarity_with": language2, "similar_word": "appā (அப்பா)"}
    word_pairs["Broom"] = {"language": language1, "word": "bitjaru (빗자루)", "similarity_with": language2, "similar_word": "thudaippam (துடைப்பம்)"}

    print(f"The two languages are {language1} and {language2}.")
    print("--------------------------------------------------")
    print("Here is a comparison of the words for 'mom', 'dad', and 'broom':\n")

    for concept, words in word_pairs.items():
        print(f"Word for \"{concept}\":")
        print(f"  - In {words['language']}: {words['word']}")
        print(f"  - In {words['similarity_with']}: {words['similar_word']}\n")

    print("Analysis:")
    print("The words for 'mom' (eomma/ammā) and 'dad' (appa/appā) are nearly identical.")
    print("This is a well-known linguistic coincidence. The words for 'broom' are not similar.")

find_similar_words()
