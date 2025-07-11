import operator

def solve_medical_case():
    """
    Analyzes patient data to find the most likely tick-borne illness.
    """
    # Patient clinical clues extracted from the case description
    patient_clues = {
        'location': 'Oklahoma',
        'exposure': 'Tick',
        'symptoms': ['fever', 'headache', 'myalgia'],
        'specific_symptoms': ['disorientation'],
        'labs': 'Negative Lyme IgG'
    }

    # Database of potential diagnoses and their typical features for scoring
    # Scores are based on how strongly a feature points to a diagnosis.
    # Geo: Geographic fit (Oklahoma)
    # Vec: Vector fit (Tick)
    # Symp: Core symptoms fit
    # Neuro: Neurologic symptoms fit
    # Labs: Consistency with lab results
    diagnoses = {
        'Babesia microti':    {'Geo': 1, 'Vec': 2, 'Symp': 2, 'Neuro': 0, 'Labs': 1},
        'Plasmodium':         {'Geo': 0, 'Vec': 0, 'Symp': 1, 'Neuro': 1, 'Labs': 1},
        'Borrelia burgdorferi':{'Geo': 1, 'Vec': 2, 'Symp': 2, 'Neuro': 1, 'Labs': -2}, # Penalized for negative IgG
        'Ehrlichia':          {'Geo': 3, 'Vec': 2, 'Symp': 2, 'Neuro': 2, 'Labs': 1}, # Strong fit for all criteria
        'Rickettsia rickettsii':{'Geo': 3, 'Vec': 2, 'Symp': 2, 'Neuro': 2, 'Labs': 1}  # Also a strong fit
    }

    # In a tie between Ehrlichia and Rickettsia, the classic presentation of fever,
    # headache, and neurologic changes without a mentioned rash slightly favors Ehrlichia.
    # We will add a small tie-breaker point.
    diagnoses['Ehrlichia']['Tie-breaker'] = 1
    diagnoses['Rickettsia rickettsii']['Tie-breaker'] = 0


    print("Evaluating potential diagnoses based on patient's case:\n")
    final_scores = {}

    for disease, criteria_scores in diagnoses.items():
        total_score = 0
        equation_parts = []
        for criterion, score in criteria_scores.items():
            total_score += score
            equation_parts.append(f"{criterion} ({score})")
        
        final_scores[disease] = total_score
        
        # Print the scoring equation for each disease
        equation_str = " + ".join(equation_parts)
        print(f"Scoring for {disease}:")
        print(f"{equation_str} = {total_score}\n")

    # Determine the most likely diagnosis
    best_diagnosis = max(final_scores.items(), key=operator.itemgetter(1))[0]

    print("-" * 30)
    print(f"Conclusion: The patient's presentation (fever, headache, myalgia, disorientation after a trip to Oklahoma) is most consistent with an infection by {best_diagnosis}.")
    print("This tick-borne illness is endemic to the region and matches the acute symptoms.")
    print("-" * 30)


solve_medical_case()