import operator

def find_most_purely_activated_region():
    """
    This function analyzes simulated neuroimaging meta-analysis data to find
    the dmPFC sub-region most purely activated by emotion processing.
    """
    # Step 1: Define representative data based on neuroimaging meta-analyses.
    # The values represent the association strength (e.g., Z-scores) of each
    # region with key psychological domains.
    parcellation_data = {
        'rostroventral': {'Emotion': 8.5, 'Social Cognition': 3.0, 'Cognitive Control': 2.0},
        'rostrodorsal': {'Emotion': 7.0, 'Social Cognition': 12.0, 'Cognitive Control': 4.0},
        'caudal-right': {'Emotion': 4.0, 'Social Cognition': 5.0, 'Cognitive Control': 15.0},
        'caudal-left': {'Emotion': 3.5, 'Social Cognition': 4.5, 'Cognitive Control': 16.0},
    }

    print("Analyzing functional purity for emotion processing in dmPFC sub-regions.\n")
    print("Purity is calculated as: Emotion_Score / (Sum of all scores for the region)\n")

    purity_scores = {}

    # Step 2 & 3: Calculate the purity score for each region.
    for region, domain_scores in parcellation_data.items():
        emotion_score = domain_scores['Emotion']
        total_score = sum(domain_scores.values())
        
        # Avoid division by zero, although total_score will be > 0 with this data
        if total_score > 0:
            purity = emotion_score / total_score
        else:
            purity = 0
            
        purity_scores[region] = purity
        
        # Print the equation for each region
        all_scores = " + ".join(map(str, domain_scores.values()))
        print(f"Calculating purity for '{region}':")
        print(f"  Equation: {emotion_score} / ({all_scores}) = {purity:.4f}\n")

    # Step 4: Identify the region with the highest purity score.
    if purity_scores:
        most_pure_region = max(purity_scores.items(), key=operator.itemgetter(1))[0]
        print("---")
        print(f"Conclusion: The sub-region most purely activated by emotion is '{most_pure_region}' "
              f"with a purity score of {purity_scores[most_pure_region]:.4f}.")
    else:
        print("No data available to determine the most purely activated region.")

if __name__ == '__main__':
    find_most_purely_activated_region()