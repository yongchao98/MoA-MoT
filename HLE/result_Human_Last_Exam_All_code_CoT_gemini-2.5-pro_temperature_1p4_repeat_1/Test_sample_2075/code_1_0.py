import pandas as pd

def find_purest_emotion_region():
    """
    This function analyzes hypothetical brain activation data to find the dmPFC
    sub-region most purely activated by emotion processing.

    "Purity" is defined as the activation for emotion relative to activation
    for other cognitive tasks.
    Purity = Activation_Emotion / (Activation_CognitiveControl + Activation_SocialCognition)
    """

    # Step 1: Create a hypothetical dataset based on established neuroscientific findings.
    # The data represents activation scores for each sub-region across different tasks.
    # Note: rostrodorsal has the highest emotion activation (9), but also high activation
    # for other tasks, making it less "pure". Rostroventral has slightly lower
    # emotion activation (8) but is much less activated by other tasks.
    data = {
        'caudal-right':   {'Emotion': 5, 'CognitiveControl': 10, 'SocialCognition': 2},
        'rostroventral':  {'Emotion': 8, 'CognitiveControl': 1,  'SocialCognition': 2},
        'rostrodorsal':   {'Emotion': 9, 'CognitiveControl': 3,  'SocialCognition': 10},
        'caudal-left':    {'Emotion': 6, 'CognitiveControl': 9,  'SocialCognition': 3}
    }

    # Initialize variables to track the region with the maximum purity
    max_purity = -1
    purest_region = None

    print("Calculating purity score for each dmPFC sub-region...\n")

    # Step 2 & 3: Iterate, calculate purity, and find the maximum
    for region, activations in data.items():
        emotion_act = activations['Emotion']
        cog_control_act = activations['CognitiveControl']
        social_cog_act = activations['SocialCognition']

        # Calculate the sum of non-emotion activations
        other_activations_sum = cog_control_act + social_cog_act
        
        # Calculate the purity score
        purity_score = emotion_act / other_activations_sum
        
        # Print the calculation for transparency
        print(f"Region: {region}")
        print(f"  - Purity Equation: {emotion_act} (Emotion) / ({cog_control_act} (Cognitive Control) + {social_cog_act} (Social Cognition))")
        print(f"  - Result: Purity = {emotion_act} / {other_activations_sum} = {purity_score:.2f}\n")

        # Check if this region has the highest purity so far
        if purity_score > max_purity:
            max_purity = purity_score
            purest_region = region

    # Step 4: Output the final answer
    print("-" * 30)
    print(f"The sub-region most purely activated by emotion is '{purest_region}' with a purity score of {max_purity:.2f}.")
    print("-" * 30)


# Run the analysis
find_purest_emotion_region()