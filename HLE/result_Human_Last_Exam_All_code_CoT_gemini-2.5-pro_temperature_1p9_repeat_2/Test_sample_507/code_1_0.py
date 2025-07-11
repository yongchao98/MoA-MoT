import pandas as pd

def find_similar_words():
    """
    Identifies and prints a comparison of words from two Asian languages
    with no direct historical contact that show surprising similarities.
    """
    print("Finding two languages in Asia from different cultures with no direct contact,")
    print("where the words for 'mom', 'dad', and 'broom' have similarities.\n")
    print("Selected Languages: Korean (Koreanic Family) and Tamil (Dravidian Family).")
    print("These languages are linguistically unrelated and geographically distant.\n")

    # Data for comparison
    data = {
        'English': ["Mom", "Dad", "Broom"],
        'Korean': ["eomma (엄마)", "appa (아빠)", "bitjaru (빗자루)"],
        'Tamil': ["amma (அம்மா)", "appa (அப்பா)", "thudaippam (துடைப்பம்)"]
    }

    # Create a DataFrame for a clean, aligned printout.
    df = pd.DataFrame(data)

    print("Word Comparison:")
    # Print the DataFrame without the index for a cleaner look.
    print(df.to_string(index=False))

    print("\n--- Analysis ---")
    print(f"Similarity for 'Mom': Korean '{data['Korean'][0]}' and Tamil '{data['Tamil'][0]}' are nearly identical.")
    print(f"Similarity for 'Dad': Korean '{data['Korean'][1]}' and Tamil '{data['Tamil'][1]}' are identical.")
    print(f"Similarity for 'Broom': Korean '{data['Korean'][2]}' and Tamil '{data['Tamil'][2]}' are not similar.")
    print("\nConclusion: While a perfect match for a culturally specific word like 'broom' is not found,")
    print("the identical nature of the universal family words 'mom' and 'dad' in these two")
    print("unrelated languages is a remarkable linguistic coincidence.")

# Execute the function to display the results.
find_similar_words()
