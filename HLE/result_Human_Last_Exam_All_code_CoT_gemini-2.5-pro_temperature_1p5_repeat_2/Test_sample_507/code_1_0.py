import pandas as pd

def find_similar_words():
    """
    This function presents a comparison of words for "mom," "dad," and "broom"
    in Korean and Persian (Farsi), two Asian languages with no direct historical contact.
    """

    data = {
        'Concept': ['Mom', 'Dad', 'Broom'],
        'Korean': ['엄마 (eomma)', '아빠 (appa)', '빗자루 (bitjaru)'],
        'Persian (Farsi)': ['مامان (māmān)', 'بابا (bābā)', 'جارو (jâru)']
    }
    
    df = pd.DataFrame(data)

    print("Two languages from distant Asian cultures with no direct contact that share highly similar words are Korean and Persian (Farsi).\n")
    print("Here is a comparison of the words for 'mom', 'dad', and 'broom':\n")
    
    # Print the DataFrame in a formatted way
    print(df.to_string(index=False))

    print("\n\n--- Analysis of Similarities ---\n")
    print("1. Mom: Korean 'eomma' and Persian 'māmān' are both common, informal words for mother that utilize the 'm' sound typical in nursery words.")
    print("2. Dad: Korean 'appa' and Persian 'bābā' are nearly identical, both using bilabial stop sounds (/p/ and /b/) common in words for father.")
    print("3. Broom: The resemblance between Korean 'bitjaru' and Persian 'jâru' is a remarkable coincidence for a non-nursery word, with the core 'jaru' sound being almost the same.")
    
find_similar_words()