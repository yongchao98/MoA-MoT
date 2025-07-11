import sys

def find_most_purely_emotional_region():
    """
    This function models the functional specialization of dmPFC sub-regions
    to find the one most purely activated by emotion processing.
    
    The scores are illustrative and based on established neuroscientific findings where:
    - rostroventral mPFC is primarily for emotion/affect.
    - rostrodorsal mPFC is for cognitive social processes (e.g., mentalizing).
    - caudal mPFC is for cognitive control and action monitoring.
    """
    # Database of functional profiles. 'Emotion' refers to affective processing,
    # 'Cognition' refers to control, mentalizing, and other non-affective processes.
    functional_profiles = {
        'caudal-right': {'Emotion': 15, 'Cognition': 85},
        'rostroventral': {'Emotion': 80, 'Cognition': 20},
        'rostrodorsal': {'Emotion': 40, 'Cognition': 60},
        'caudal-left': {'Emotion': 15, 'Cognition': 85},
    }

    # To illustrate the "strength vs. purity" point made in the question:
    # A region could have stronger emotion activation (e.g., 90) but be less pure
    # if it's also highly active for cognition.
    # For example: {'Emotion': 90, 'Cognition': 50} -> Purity = 90 / (90+50) = 0.64
    # This is less pure than the 'rostroventral' region's score of 0.80.

    most_pure_region = None
    highest_purity_score = -1.0

    print("Calculating 'Emotion Purity' score for each dmPFC sub-region.")
    print("Purity Score = Emotion Score / (Emotion Score + Cognition Score)\n")

    for region, scores in functional_profiles.items():
        emotion_score = scores['Emotion']
        cognition_score = scores['Cognition']
        
        # Avoid division by zero, although not expected with this data
        if (emotion_score + cognition_score) == 0:
            purity_score = 0.0
        else:
            # Calculate the purity score
            purity_score = emotion_score / (emotion_score + cognition_score)

        print(f"Region: {region}")
        # Here we output each number in the final equation
        print(f"  Calculation: {emotion_score} / ({emotion_score} + {cognition_score}) = {purity_score:.2f}")

        # Check if this region is the new "most pure"
        if purity_score > highest_purity_score:
            highest_purity_score = purity_score
            most_pure_region = region
        print("-" * 30)

    print(f"\nResult:")
    print(f"The sub-region most purely activated by emotion processing is '{most_pure_region}'.")

# Execute the function
find_most_purely_emotional_region()