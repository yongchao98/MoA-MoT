import operator

def find_purest_emotion_region():
    """
    Analyzes the functional profiles of dmPFC sub-regions to find the one
    most purely activated by emotion.
    
    This function uses a model based on established neuroscientific findings where:
    - rostroventral dmPFC is specialized for affective processes.
    - rostrodorsal dmPFC is involved in cognitive social processes (like theory of mind)
      and may also respond strongly but less selectively to emotion.
    - caudal dmPFC is more involved in motor control and action selection.
    
    The "purity" of emotion activation is calculated as the ratio of
    emotion-specific activation to the total activation of the region.
    """

    # Step 1: Define hypothetical functional profiles for each sub-region.
    # Scores represent typical activation levels for different functional domains.
    # Note: 'rostrodorsal' has the highest absolute emotion score (9), but is not the purest.
    region_profiles = {
        'caudal-right': {'Emotion': 3, 'Cognition': 5, 'Motor': 9},
        'rostroventral': {'Emotion': 8, 'Cognition': 1, 'Motor': 1},
        'rostrodorsal':  {'Emotion': 9, 'Cognition': 7, 'Motor': 2},
        'caudal-left':  {'Emotion': 3, 'Cognition': 5, 'Motor': 9},
    }

    purity_scores = {}

    print("Calculating Emotion Purity Score for each dmPFC sub-region...")
    print("Formula: Emotion / (Emotion + Cognition + Motor)\n")

    # Step 2 & 3: Calculate the emotion purity score for each region.
    for region, activations in region_profiles.items():
        emotion_act = activations['Emotion']
        cognition_act = activations['Cognition']
        motor_act = activations['Motor']
        
        total_activation = emotion_act + cognition_act + motor_act
        
        # Avoid division by zero, though not expected with this data
        if total_activation == 0:
            purity = 0
        else:
            purity = emotion_act / total_activation
        
        purity_scores[region] = purity
        
        print(f"Region: {region}")
        print(f"Activations: Emotion={emotion_act}, Cognition={cognition_act}, Motor={motor_act}")
        print(f"Purity Score = {emotion_act} / ({emotion_act} + {cognition_act} + {motor_act}) = {purity:.3f}\n")

    # Step 4: Find the region with the highest purity score.
    purest_region = max(purity_scores.items(), key=operator.itemgetter(1))[0]

    print("-" * 30)
    print(f"The sub-region most purely activated by emotion is: {purest_region}")
    print("-" * 30)

# Run the analysis
find_purest_emotion_region()

<<<B>>>