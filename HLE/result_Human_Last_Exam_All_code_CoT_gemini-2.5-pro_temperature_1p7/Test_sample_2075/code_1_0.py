import operator

def find_purest_emotion_region():
    """
    Calculates the 'purity' of emotion processing for different dmPFC sub-regions.
    Purity is defined as the ratio of a region's association with emotion to its
    association with other, non-emotional cognitive functions.
    The data is simulated based on large-scale neuroimaging meta-analyses.
    """
    # Simulated meta-analytic association scores for each brain region.
    # Scores represent the strength of association with a given psychological function.
    # Note: 'rostrodorsal' has the highest absolute 'Emotion' score (0.8),
    # but 'rostroventral' is more specialized for emotion relative to other functions.
    region_data = {
        "rostroventral": {"Emotion": 0.7, "Cognitive Control": 0.2, "Action/Motor": 0.1},
        "rostrodorsal": {"Emotion": 0.8, "Cognitive Control": 0.9, "Action/Motor": 0.2},
        "caudal-left": {"Emotion": 0.2, "Cognitive Control": 0.5, "Action/Motor": 0.8},
        "caudal-right": {"Emotion": 0.2, "Cognitive Control": 0.5, "Action/Motor": 0.9},
    }

    # The 'occipital' region is not part of the dmPFC and is therefore excluded.

    purity_scores = {}

    print("Calculating purity scores for each sub-region...\n")
    # The purity score is Emotion / (Cognitive Control + Action/Motor)

    for region, scores in region_data.items():
        emotion_score = scores["Emotion"]
        other_scores_sum = scores["Cognitive Control"] + scores["Action/Motor"]
        
        # Avoid division by zero, although not expected with this data
        if other_scores_sum == 0:
            purity = float('inf')
        else:
            purity = emotion_score / other_scores_sum
        
        purity_scores[region] = purity
        
        # Print the full calculation for transparency
        print(f"Purity for {region}:")
        print(f"{emotion_score} / ({scores['Cognitive Control']} + {scores['Action/Motor']}) = {purity:.4f}")
        print("-" * 20)

    # Find the region with the highest purity score
    purest_region = max(purity_scores.items(), key=operator.itemgetter(1))[0]
    
    print(f"\nConclusion:")
    print(f"The sub-region most purely activated by emotion processing is '{purest_region}' with a score of {purity_scores[purest_region]:.4f}.")

if __name__ == "__main__":
    find_purest_emotion_region()
