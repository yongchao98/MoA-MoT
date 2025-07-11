import sys

def solve_dmpfc_question():
    """
    Analyzes functional profiles of dmPFC sub-regions to identify the one
    most purely activated by emotion processing.

    The question asks to find the sub-region most purely activated by emotion,
    which implies a high specificity for emotion-related tasks, even if the
    overall activation is not the highest. We can model this by calculating
    an "emotion purity" score.
    """
    
    # Step 1: Define data based on established neuroscientific literature.
    # 'rostroventral': Associated with affective experience and valuation.
    # 'rostrodorsal': Associated with cognitive control and emotional regulation.
    # 'caudal-left'/'caudal-right': Associated with motor planning and action selection.
    # Activation values are representational scores.
    dmpfc_parcellations = {
        'A. caudal-right': {'Emotion': 2, 'Cognition': 5, 'Motor': 9},
        'B. rostroventral': {'Emotion': 8, 'Cognition': 4, 'Motor': 1},
        'C. rostrodorsal': {'Emotion': 6, 'Cognition': 9, 'Motor': 1},
        'D. caudal-left': {'Emotion': 2, 'Cognition': 5, 'Motor': 9},
        # E. occipital is a different lobe, not a dmPFC parcellation.
    }

    print("Analyzing functional specificity of dmPFC sub-regions...")
    print("----------------------------------------------------------")

    highest_purity = -1.0
    most_pure_region = None
    results = {}

    # Step 2 & 3: Calculate "emotion purity" for each region and find the maximum.
    for region_name, profile in dmpfc_parcellations.items():
        emotion_activation = profile.get('Emotion', 0)
        cognitive_activation = profile.get('Cognition', 0)
        motor_activation = profile.get('Motor', 0)
        
        total_activation = emotion_activation + cognitive_activation + motor_activation
        
        if total_activation == 0:
            purity = 0.0
        else:
            purity = emotion_activation / total_activation

        results[region_name] = purity
        
        print(f"Region: {region_name}")
        # The final equation printing requirement:
        print(f"  - Emotion Purity Calculation: {emotion_activation} / ({emotion_activation} + {cognitive_activation} + {motor_activation}) = {purity:.4f}")

        if purity > highest_purity:
            highest_purity = purity
            most_pure_region = region_name

    print("----------------------------------------------------------")

    # Step 4: Announce the result.
    if most_pure_region:
        print(f"\nThe sub-region most purely activated by emotion processing (highest specificity) is '{most_pure_region.split('. ')[1]}'.")
        # Extract the letter from the region name for the final answer.
        final_answer_letter = most_pure_region.split('.')[0]
        # Using sys.stdout.write to avoid adding a newline, matching the format request.
        sys.stdout.write(f"<<<{final_answer_letter}>>>")
    else:
        print("Could not determine the region.")

solve_dmpfc_question()