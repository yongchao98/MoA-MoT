def find_linguistic_coincidences():
    """
    This function identifies and displays a comparison between Korean and Telugu,
    two languages from different Asian cultures with no direct contact, which
    share surprisingly similar words for "mom", "dad", and "broom".
    """
    
    print("Finding two languages from different Asian cultures with no direct contact")
    print("that have similar words for 'mom', 'dad', and 'broom'.\n")

    # Define the languages and their details
    language1_name = "Korean"
    language2_name = "Telugu"
    
    # Store the word data in a dictionary for easy access and printing
    comparison_data = {
        "Word": ("Mom", "Dad", "Broom"),
        "Korean": ("eomma (엄마)", "appa (아빠)", "bitjaru (빗자루)"),
        "Telugu": ("amma (అమ్మ)", "appa (అప్ప)", "chīpuru (చీపురు)")
    }
    
    # Print the findings
    print(f"Language 1: {language1_name} (Koreanic Family, East Asia)")
    print(f"Language 2: {language2_name} (Dravidian Family, South Asia)")
    print("\nThese two languages are unrelated and their cultures developed without direct contact.\n")
    
    print("--- Word Comparison ---")
    header = f"{'Word':<10} | {language1_name:<20} | {language2_name:<20}"
    print(header)
    print("-" * len(header))
    
    # Loop through the words and print the comparison
    for i in range(len(comparison_data["Word"])):
        word = comparison_data["Word"][i]
        word_lang1 = comparison_data[language1_name][i]
        word_lang2 = comparison_data[language2_name][i]
        print(f"{word:<10} | {word_lang1:<20} | {word_lang2:<20}")
        
    print("-" * len(header))
    
    print("\nAnalysis of Similarity:")
    print("- Mom (eomma/amma) and Dad (appa/appa) are nearly identical.")
    print("- Broom (bitjaru/chīpuru) shows coincidental structural similarity,")
    print("  especially with the 'i' sound in the first syllable and the '-ru' ending.")

# Execute the function to print the result
find_linguistic_coincidences()