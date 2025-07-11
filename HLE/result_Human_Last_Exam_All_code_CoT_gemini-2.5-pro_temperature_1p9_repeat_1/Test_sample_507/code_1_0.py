import pandas as pd

def find_similar_words():
    """
    Identifies and displays similar words in Korean and Tamil, two languages
    from disparate Asian cultures.
    """
    # Data for the comparison
    data = {
        'Term': ["Mom", "Dad", "Broom"],
        'Korean': ["Eomma (엄마)", "Appa (아빠)", "Bitjaru (빗자루)"],
        'Tamil': ["Amma (அம்மா)", "Appa (அப்பா)", "Thudaippam (துடைப்பம்)"]
    }

    # Create a DataFrame for clean, aligned printing
    df = pd.DataFrame(data)

    # Print the introductory statement and the DataFrame
    print("The two languages are Korean and Tamil.")
    print("These represent distinct cultures in East and South Asia, respectively, with no direct contact.")
    print("-" * 65)
    print("Here is a comparison of the words for 'mom', 'dad', and 'broom':\n")
    print(df.to_string(index=False))
    print("-" * 65)
    print("\nAs shown, the words for 'Mom' and 'Dad' are nearly identical,")
    print("a remarkable similarity for two languages from different families.")

if __name__ == '__main__':
    find_similar_words()