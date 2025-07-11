def find_correct_linguistic_features():
    """
    Analyzes choices to find the correct group of linguistic features for word complexity.
    """
    choices = {
        "A": "word length, word familiarity rating, number of unique words, concreteness rating",
        "B": "word familiarity rating, number of word categories, syllable count, concreteness rating",
        "C": "word familiarity rating, syllable count, concreteness rating, number of unique words",
        "D": "word length, word familiarity rating, syllable count, concreteness rating",
        "E": "word length, imageability rating, word familiarity rating, number of word categories"
    }

    # The question is about the complexity of individual *words*.
    # Therefore, features that describe a *text* are not applicable.
    text_level_features = ["number of unique words", "number of word categories"]

    print("Analyzing the choices to find the one with only word-level linguistic features...")
    print("-" * 70)

    correct_choice = None
    for choice, features_text in choices.items():
        contains_text_level_feature = False
        for feature in text_level_features:
            if feature in features_text:
                print(f"Choice {choice} is eliminated because it contains the text-level feature: '{feature}'.")
                contains_text_level_feature = True
                break
        
        if not contains_text_level_feature:
            correct_choice = choice
            print(f"Choice {choice} is the correct answer. It contains only word-level features: '{features_text}'.")

    print("-" * 70)
    
    if correct_choice:
        # The prompt asks to output the final answer in a specific format.
        # The final answer is the letter of the correct choice.
        print(f"\n<<<D>>>")

find_correct_linguistic_features()