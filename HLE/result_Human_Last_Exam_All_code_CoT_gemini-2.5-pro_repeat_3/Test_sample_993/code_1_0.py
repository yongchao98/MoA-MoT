def analyze_paleoclimate_factors():
    """
    Analyzes environmental factors to determine which conditions lead to low
    levels of the 18O isotope in lakebed sulfate deposits.
    """
    print("Step 1: Define the scientific principle.")
    print("Low levels of the 18O isotope are caused by low rates of evaporation.")
    print("We will create a scoring model where factors that decrease evaporation get a negative score.\n")

    # A score of -1 means the factor contributes to LOW 18O levels.
    # A score of +1 means the factor contributes to HIGH 18O levels.
    factor_scores = {
        'Climate': {'Wet': -1, 'Dry': 1},
        'Temperature': {'Cold': -1, 'Warm': 1},
        'Lake Level': {'High': -1, 'Shallow': 1}
    }

    # The provided answer choices
    choices = {
        'A': {'Climate': 'Wet', 'Temperature': 'Warm', 'Lake Level': 'Shallow'},
        'B': {'Climate': 'Dry', 'Temperature': 'Warm', 'Lake Level': 'Shallow'},
        'C': {'Climate': 'Wet', 'Temperature': 'Cold', 'Lake Level': 'Shallow'},
        'D': {'Climate': 'Dry', 'Temperature': 'Cold', 'Lake Level': 'Shallow'},
        'E': {'Climate': 'Wet', 'Temperature': 'Warm', 'Lake Level': 'High'},
        'F': {'Climate': 'Dry', 'Temperature': 'Warm', 'Lake Level': 'High'},
        'G': {'Climate': 'Wet', 'Temperature': 'Cold', 'Lake Level': 'High'},
        'H': {'Climate': 'Dry', 'Temperature': 'Cold', 'Lake Level': 'High'}
    }

    print("Step 2: Calculate an '18O Enrichment Score' for each choice.")
    print("The score is the sum of the factor scores. The lowest score is the answer.\n")
    
    results = {}
    for choice, factors in choices.items():
        climate = factors['Climate']
        temp = factors['Temperature']
        level = factors['Lake Level']

        # Get the score for each factor
        climate_score = factor_scores['Climate'][climate]
        temp_score = factor_scores['Temperature'][temp]
        level_score = factor_scores['Lake Level'][level]

        # Calculate the total score
        total_score = climate_score + temp_score + level_score
        results[choice] = total_score
        
        # Print the equation for each choice
        print(f"Choice {choice}: {climate}({climate_score}) + {temp}({temp_score}) + {level}({level_score}) = {total_score}")

    # Find the choice with the minimum score
    best_choice = min(results, key=results.get)

    print("\nStep 3: Conclusion.")
    print(f"The analysis shows that choice '{best_choice}' has the lowest score, indicating it represents the conditions most responsible for low levels of the 18O isotope.")
    print("These conditions are: a Wet, Cold climate with High lake levels.")

if __name__ == '__main__':
    analyze_paleoclimate_factors()