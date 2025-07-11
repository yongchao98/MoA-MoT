import pandas as pd

def analyze_linguistic_features():
    """
    Analyzes and identifies the most relevant group of linguistic features for
    predicting word complexity based on established psycholinguistic principles.
    """

    # Principle 1: Word-level features are more predictive of word complexity
    # than text-level features. We assign higher scores to word-level features.
    # Scores are heuristic, representing importance.
    feature_relevance = {
        # Core word-level predictors
        'word familiarity rating': 3,  # Strongest predictor (frequency/exposure)
        'syllable count': 3,         # Strong predictor (phonological/structural complexity)
        'concreteness rating': 3,    # Strong predictor (semantic processing ease)
        'word length': 2,            # Good predictor, but often correlated with syllable count
        'imageability rating': 2,    # Good predictor, but highly correlated with concreteness

        # Features more related to text-level complexity
        'number of unique words': 1,  # Text-level (lexical diversity)
        'number of word categories': 1 # Text-level (syntactic variety)
    }

    # The provided answer choices
    choices = {
        'A': ['word length', 'word familiarity rating', 'number of unique words', 'concreteness rating'],
        'B': ['word familiarity rating', 'number of word categories', 'syllable count', 'concreteness rating'],
        'C': ['word familiarity rating', 'syllable count', 'concreteness rating', 'number of unique words'],
        'D': ['word length', 'word familiarity rating', 'syllable count', 'concreteness rating'],
        'E': ['word length', 'imageability rating', 'word familiarity rating', 'number of word categories']
    }

    results = []
    print("Evaluating feature combinations for predicting word complexity:")
    print("-" * 60)

    for choice, features in choices.items():
        score = 0
        equation_parts = []
        for feature in features:
            relevance_score = feature_relevance.get(feature, 0)
            score += relevance_score
            equation_parts.append(str(relevance_score))
        
        # This fulfills the "output each number in the final equation" requirement
        equation_str = " + ".join(equation_parts)
        
        results.append({
            'Choice': choice,
            'Features': ", ".join(features),
            'Score Equation': f"{equation_str} = {score}",
            'Score': score
        })

    # Display results in a structured way
    results_df = pd.DataFrame(results)
    print(results_df[['Choice', 'Score Equation', 'Features']].to_string(index=False))
    print("-" * 60)

    # Determine the best choice
    best_choice_row = results_df.loc[results_df['Score'].idxmax()]
    best_choice = best_choice_row['Choice']

    print(f"Conclusion: Choice {best_choice} has the highest score.")
    print("It combines the most powerful and distinct word-level predictors: measures of structure (length/syllables), frequency (familiarity), and meaning (concreteness).")

analyze_linguistic_features()
<<<D>>>