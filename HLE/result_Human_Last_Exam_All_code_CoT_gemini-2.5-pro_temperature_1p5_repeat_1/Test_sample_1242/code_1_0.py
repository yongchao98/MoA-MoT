def calculate_word_complexity(word_data):
    """
    Calculates a hypothetical complexity score for a word based on its linguistic features.
    
    Args:
        word_data (dict): A dictionary containing the word and its features.
    """
    
    # These are the key features from the correct answer choice
    length = word_data['length']
    syllables = word_data['syllables']
    familiarity = word_data['familiarity_rating'] # Assumed on a 1-7 scale (higher is more familiar)
    concreteness = word_data['concreteness_rating'] # Assumed on a 1-7 scale (higher is more concrete)
    
    # Define weights for each feature.
    # Positive weights for features that increase complexity (length, syllables).
    # Negative weights for features that decrease complexity (familiarity, concreteness).
    w_len = 0.1
    w_syl = 0.5
    w_fam = -0.4
    w_con = -0.3

    # Calculate the complexity score
    score = (w_len * length) + (w_syl * syllables) + (w_fam * familiarity) + (w_con * concreteness)
    
    # Print the equation with the actual numbers, as requested
    print(f"Calculating complexity for the word: '{word_data['word']}'")
    print("Features: word length, syllable count, familiarity rating, concreteness rating.")
    print("\nHypothetical Complexity Equation:")
    print(f"({w_len} * {length}) + ({w_syl} * {syllables}) + ({w_fam} * {familiarity}) + ({w_con} * {concreteness}) = {score:.2f}")

# --- Main execution ---
# Example data for a complex word
complex_word = {
    'word': 'ambiguous',
    'length': 9,
    'syllables': 4,
    'familiarity_rating': 3.5, # Moderately unfamiliar
    'concreteness_rating': 2.1  # Highly abstract
}

# Example data for a simple word
simple_word = {
    'word': 'apple',
    'length': 5,
    'syllables': 2,
    'familiarity_rating': 6.8, # Highly familiar
    'concreteness_rating': 6.9  # Highly concrete
}

# Run the calculation and print the results for both words
calculate_word_complexity(complex_word)
print("-" * 30)
calculate_word_complexity(simple_word)
