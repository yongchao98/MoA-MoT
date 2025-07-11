import pandas as pd

def find_most_selective_region():
    """
    This function simulates an analysis to find the dmPFC sub-region most
    selectively activated by emotion processing.

    The 'purity' of activation is defined as the ratio of emotion activation
    to the sum of activations for other major cognitive functions.
    A higher score indicates greater selectivity for emotion.
    """
    # Step 1: Define hypothetical activation data based on neuroscience literature.
    # Scores represent the strength of association of a region with a function.
    data = {
        'region': ['caudal-right', 'rostroventral', 'rostrodorsal', 'caudal-left', 'occipital'],
        'emotion_activation': [3.0, 9.0, 4.0, 3.0, 1.0],  # rostroventral is strongly tied to emotion/valuation
        'cognition_activation': [5.0, 2.0, 9.0, 5.0, 1.0], # rostrodorsal is for theory-of-mind/metacognition
        'control_activation': [8.0, 1.0, 2.0, 9.0, 1.0]   # caudal regions are for cognitive/motor control
    }
    df = pd.DataFrame(data)

    # Step 2: Calculate the "purity" score for emotion activation.
    # Purity = Emotion / (Cognition + Control)
    # We add a small epsilon (0.01) to the denominator to avoid division by zero.
    df['other_activation_sum'] = df['cognition_activation'] + df['control_activation']
    df['emotion_purity_score'] = df['emotion_activation'] / (df['other_activation_sum'] + 0.01)

    # Step 3: Find the region with the maximum purity score.
    most_selective_region = df.loc[df['emotion_purity_score'].idxmax()]

    # Step 4: Print the results and the final calculation.
    print("Analysis of dmPFC Sub-region Selectivity for Emotion Processing\n")
    print(df[['region', 'emotion_purity_score']].round(2))
    print("\n--- Finding the Most Purely Activated Region ---")

    winner_name = most_selective_region['region']
    emotion_score = most_selective_region['emotion_activation']
    cognition_score = most_selective_region['cognition_activation']
    control_score = most_selective_region['control_activation']
    purity_score = most_selective_region['emotion_purity_score']

    print(f"\nThe region most purely activated by emotion is: '{winner_name}'")
    print(f"\nThis conclusion is based on the highest emotion purity score.")
    print("\nFinal Equation for the winning region:")
    print(f"Purity Score = Emotion Activation / (Cognition Activation + Control Activation)")
    print(f"'{winner_name}' Purity Score = {emotion_score} / ({cognition_score} + {control_score}) = {purity_score:.2f}")

find_most_selective_region()