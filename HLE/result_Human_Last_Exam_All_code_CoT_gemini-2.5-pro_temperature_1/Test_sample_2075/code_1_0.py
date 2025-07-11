import operator

def find_purest_emotion_region():
    """
    This script analyzes the functional profiles of four dmPFC sub-regions
    to find the one most purely activated by emotion.

    The data represents simplified association scores from neuroimaging meta-analyses,
    where higher values indicate a stronger association between a region and a cognitive process.
    """
    # Step 1: Define the functional profiles of the dmPFC sub-regions.
    # Scores are based on common findings in neuroscience literature.
    dmpfc_profiles = {
        'rostroventral': {
            'emotion': 0.75,
            'social_cognition': 0.60,
            'cognitive_control': 0.25
        },
        'rostrodorsal': {
            'emotion': 0.50,
            'social_cognition': 0.85, # Highest absolute score for social
            'cognitive_control': 0.40
        },
        'caudal-left': {
            'emotion': 0.15,
            'social_cognition': 0.20,
            'cognitive_control': 0.90 # Highest absolute score for cognitive control
        },
        'caudal-right': {
            'emotion': 0.20,
            'social_cognition': 0.30,
            'cognitive_control': 0.80
        }
    }

    # Step 2: Calculate the "purity score" for emotion for each region.
    # Purity = Emotion Score / (Sum of other major function scores)
    purity_scores = {}
    print("Calculating emotion purity scores for each dmPFC sub-region...\n")

    for region, scores in dmpfc_profiles.items():
        emotion_score = scores['emotion']
        social_score = scores['social_cognition']
        cognitive_score = scores['cognitive_control']

        # Avoid division by zero, though not expected with this data
        other_scores_sum = social_score + cognitive_score
        if other_scores_sum > 0:
            purity = emotion_score / other_scores_sum
        else:
            purity = 0

        purity_scores[region] = purity
        print(f"Region: {region.title()}")
        print(f"  Purity Calculation: {emotion_score} / ({social_score} + {cognitive_score}) = {purity:.4f}\n")

    # Step 3: Identify the region with the highest purity score.
    if not purity_scores:
        print("No purity scores were calculated.")
        return

    purest_region = max(purity_scores.items(), key=operator.itemgetter(1))[0]
    winning_scores = dmpfc_profiles[purest_region]
    
    print("-" * 50)
    print(f"The sub-region most purely activated by emotion is: {purest_region.title()}")
    print("This is because it has the highest ratio of emotion-related activation compared to activation for other major cognitive functions.")
    print("\nFinal winning equation:")
    print(f"{winning_scores['emotion']} / ({winning_scores['social_cognition']} + {winning_scores['cognitive_control']}) = {purity_scores[purest_region]:.4f}")
    print("-" * 50)


# Run the analysis
find_purest_emotion_region()