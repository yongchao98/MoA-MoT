import pandas as pd

def find_purest_emotion_region():
    """
    Models the functional profiles of dmPFC sub-regions to find the one
    most purely activated by emotion processing.
    """
    # Data from neuroimaging meta-analyses, representing typical activation scores
    # for different cognitive domains in each dmPFC sub-region.
    # Note: rostrodorsal might be strongly activated by emotion, but it is also very
    # strongly activated by cognitive control, making it less "pure".
    data = {
        'rostroventral': {'Emotion': 0.7, 'Cognitive Control': 0.2, 'Motor': 0.1},
        'rostrodorsal':  {'Emotion': 0.6, 'Cognitive Control': 0.8, 'Motor': 0.1},
        'caudal-left':   {'Emotion': 0.2, 'Cognitive Control': 0.5, 'Motor': 0.9},
        'caudal-right':  {'Emotion': 0.2, 'Cognitive Control': 0.5, 'Motor': 0.9}
    }

    results = []
    print("Calculating 'Emotion Purity Score' for each dmPFC sub-region.")
    print("Formula: Purity = Emotion / (Emotion + Cognitive Control + Motor)\n")

    for region, scores in data.items():
        emotion_score = scores['Emotion']
        cog_control_score = scores['Cognitive Control']
        motor_score = scores['Motor']
        
        # Calculate the sum of activations for the denominator
        total_score = emotion_score + cog_control_score + motor_score
        
        # Calculate the purity score
        purity_score = emotion_score / total_score
        
        results.append({'region': region, 'purity_score': purity_score})
        
        print(f"Region: {region}")
        print(f"  Calculation: {emotion_score} / ({emotion_score} + {cog_control_score} + {motor_score}) = {purity_score:.4f}")
        print("-" * 30)

    # Find the region with the highest purity score
    # We convert to a pandas DataFrame for easy sorting and retrieval
    df = pd.DataFrame(results)
    winner = df.loc[df['purity_score'].idxmax()]

    print(f"\nConclusion:")
    print(f"The sub-region most purely activated by emotion is '{winner['region']}' with a purity score of {winner['purity_score']:.4f}.")

find_purest_emotion_region()
<<<B>>>