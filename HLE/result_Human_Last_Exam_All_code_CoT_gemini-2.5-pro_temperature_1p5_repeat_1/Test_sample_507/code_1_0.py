import pandas as pd

def find_similar_language_words():
    """
    Identifies and displays similar words from two Asian languages
    that have had no significant direct modern contact.
    """
    print("Finding two languages from two separate Asian cultures with similar words for 'mom', 'dad', and 'broom'.")
    print("The chosen languages are Turkish and Mongolian.")
    print("\nJustification:")
    print("1. Cultures: Turkish culture is rooted in Anatolia (West Asia), while Mongolian culture is in East/Central Asia.")
    print("2. Contact: The peoples have been geographically and culturally separate for many centuries with no major direct contact.")
    print("3. Linguistic Similarity: Despite the separation, they share notable similarities in core vocabulary, potentially from a distant common ancestor.\n")

    # Data for the words
    # For 'mom' and 'dad', we use ancestral forms ('Ana'/'Ata') for Turkish to show the root similarity.
    data = {
        'Concept': ['Mom', 'Dad', 'Broom'],
        'Turkish': ['Ana', 'Ata', 'Süpürge'],
        'Mongolian': ['Ene', 'Aav', 'Шүүр (Shüür)'],
        'Notes': [
            "'Ana' (Turkish) and 'Ene' (Mongolian) are phonetically similar.",
            "'Ata' (ancestral Turkic) and 'Aav' (Mongolian) are phonetically similar.",
            "The root sounds 'Süpür-' and 'Shüür' share a sibilant 's/sh', a 'ü' vowel, and an 'r'."
        ]
    }
    
    df = pd.DataFrame(data)

    print("Word Comparison:")
    # Print the dataframe in a readable format without the index
    print(df.to_string(index=False))

find_similar_language_words()
