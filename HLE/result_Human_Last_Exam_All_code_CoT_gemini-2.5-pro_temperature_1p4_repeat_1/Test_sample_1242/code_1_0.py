def calculate_complexity(word, length, syllables, familiarity, concreteness):
    """
    Calculates a hypothetical complexity score for a word based on four key features.
    
    - Complexity increases with length and number of syllables.
    - Complexity decreases with familiarity and concreteness (i.e., abstract words are more complex).
    
    Args:
        word (str): The word being evaluated.
        length (int): Number of letters in the word.
        syllables (int): Number of syllables in the word.
        familiarity (float): A rating of how familiar the word is (e.g., on a 1-7 scale).
        concreteness (float): A rating of how concrete the word is (e.g., on a 1-7 scale).
    
    Returns:
        float: A calculated complexity score.
    """
    # A simple weighted formula. Higher score means more complex.
    # We add length and syllables, and subtract familiarity and concreteness.
    # Weights are added to give more importance to certain features.
    # Here we use simple weights for demonstration.
    complexity_score = (length * 0.5) + (syllables * 1.5) - (familiarity * 1.0) - (concreteness * 1.0)
    
    print(f"Calculating complexity for the word: '{word}'")
    print(f"Equation: (length * 0.5) + (syllables * 1.5) - (familiarity * 1.0) - (concreteness * 1.0)")
    print(f"Values:   ({length} * 0.5) + ({syllables} * 1.5) - ({familiarity} * 1.0) - ({concreteness} * 1.0) = {complexity_score:.2f}")
    print("-" * 30)
    return complexity_score

# --- Data for Example Words ---
# Ratings are approximate on a 1-7 scale, where 7 is high familiarity/concreteness.

# A simple, common word
word1_data = {
    "word": "sun",
    "length": 3,
    "syllables": 1,
    "familiarity": 6.8,
    "concreteness": 6.9
}

# A complex, scientific word
word2_data = {
    "word": "photosynthesis",
    "length": 14,
    "syllables": 5,
    "familiarity": 3.5,
    "concreteness": 2.1
}

# --- Calculation ---
score1 = calculate_complexity(**word1_data)
score2 = calculate_complexity(**word2_data)

print(f"Final Score for '{word1_data['word']}': {score1:.2f}")
print(f"Final Score for '{word2_data['word']}': {score2:.2f}")
