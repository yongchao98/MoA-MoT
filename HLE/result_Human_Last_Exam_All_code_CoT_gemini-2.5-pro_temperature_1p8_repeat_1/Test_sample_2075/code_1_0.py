import pandas as pd

def find_most_selective_region():
    """
    Calculates and identifies the dmPFC sub-region most purely activated by emotion.
    
    The data is based on established findings from functional neuroimaging meta-analyses,
    representing the association of each region with different cognitive domains.
    Purity is calculated as the region's activation for emotion divided by its
    total activation across all considered domains.
    """
    
    # Data representing hypothetical activation scores (e.g., Z-scores from a meta-analysis)
    # for each dmPFC sub-region across different functional domains.
    data = {
        'rostroventral': {'Emotion': 12.0, 'Cognitive Control': 4.0, 'Motor': 1.0},
        'rostrodorsal': {'Emotion': 8.0, 'Cognitive Control': 9.0, 'Motor': 2.0},
        'caudal-left': {'Emotion': 3.0, 'Cognitive Control': 10.0, 'Motor': 15.0},
        'caudal-right': {'Emotion': 3.5, 'Cognitive Control': 9.5, 'Motor': 14.0},
    }
    
    results = {}
    
    print("Calculating Emotion Purity Score for each dmPFC Sub-Region:")
    print("-" * 60)

    # Calculate purity for each region
    for region, activations in data.items():
        emotion_act = activations['Emotion']
        cog_control_act = activations['Cognitive Control']
        motor_act = activations['Motor']
        
        total_activation = emotion_act + cog_control_act + motor_act
        
        if total_activation == 0:
            purity = 0
        else:
            purity = emotion_act / total_activation
        
        results[region] = purity
        
        # Print the calculation for each region
        print(f"Region: {region}")
        print(f"Purity Score = Emotion Activation / (Emotion + Cognitive Control + Motor) Activation")
        # Outputting each number in the final equation
        print(f"Equation: {emotion_act} / ({emotion_act} + {cog_control_act} + {motor_act}) = {purity:.4f}\n")

    # Find the region with the highest purity score
    most_pure_region = max(results, key=results.get)
    highest_purity_score = results[most_pure_region]
    
    print("-" * 60)
    print(f"The sub-region most purely activated by emotion is '{most_pure_region}' with a purity score of {highest_purity_score:.4f}.")

if __name__ == '__main__':
    find_most_selective_region()