def diagnose_ankle_pain():
    """
    Analyzes clinical findings to determine the most likely diagnosis
    by assigning scores to each possibility.
    """
    # Clinical findings and their scoring impact.
    # A positive score means the finding supports the diagnosis.
    # A large negative score means the finding strongly contradicts the diagnosis.
    findings = {
        'acute_inflammation': 2,
        'negative_xray_early': 1,
        'worsened_on_steroids': 5, # Very high impact finding
        'aspirate_no_crystals': 5,  # Rules out crystal arthropathies
        'aspirate_no_wbcs': 5      # Rules out septic/inflammatory arthritis
    }

    # Scoring matrix for each diagnosis based on the findings.
    # 1: Supports, 0: Neutral/Possible, -1: Contradicts
    diagnoses = {
        'A. Osteoarthritis': {
            'acute_inflammation': -1, # Atypical
            'negative_xray_early': 0,
            'worsened_on_steroids': -1,
            'aspirate_no_crystals': 1,
            'aspirate_no_wbcs': 1,
            'explanation': 'Presents with pain but typically not with such severe acute inflammation or worsening on steroids.'
        },
        'B. Charcot Arthropathy': {
            'acute_inflammation': 1, # Classic "hot swollen" joint
            'negative_xray_early': 1, # Classic early finding
            'worsened_on_steroids': 1, # Does not respond to anti-inflammatories
            'aspirate_no_crystals': 1,
            'aspirate_no_wbcs': 1, # Classic bland aspirate
            'explanation': 'Classic presentation: an acutely inflamed, swollen joint with bland (non-inflammatory, non-infectious) joint fluid, often triggered by minor trauma and unresponsive to steroids.'
        },
        'C. Septic Arthritis': {
            'acute_inflammation': 1,
            'negative_xray_early': 1,
            'worsened_on_steroids': 0,
            'aspirate_no_crystals': 1,
            'aspirate_no_wbcs': -1, # This finding essentially rules it out
            'explanation': 'Ruled out by the absence of white blood cells or organisms in the joint fluid.'
        },
        'D. Chronic osteomyelitis': {
            'acute_inflammation': 1,
            'negative_xray_early': -1, # Chronic process should have x-ray changes
            'worsened_on_steroids': 0,
            'aspirate_no_crystals': 1,
            'aspirate_no_wbcs': 0,
            'explanation': 'Unlikely given the negative x-rays and clean joint aspirate.'
        },
        'E. Pseudogout': {
            'acute_inflammation': 1,
            'negative_xray_early': 1,
            'worsened_on_steroids': -1, # Should improve with steroids
            'aspirate_no_crystals': -1, # This finding rules it out
            'aspirate_no_wbcs': -1, # Should have high WBCs
            'explanation': 'Ruled out by the absence of crystals in the joint fluid.'
        }
    }

    results = {}
    equations = {}
    for diagnosis, scores in diagnoses.items():
        total_score = 0
        equation_parts = []
        for finding, impact_multiplier in scores.items():
            score_value = findings[finding] * impact_multiplier
            total_score += score_value
            equation_parts.append(str(score_value))
        results[diagnosis] = total_score
        equations[diagnosis] = " + ".join(equation_parts).replace("+ -", "- ")

    # Find the best diagnosis
    best_diagnosis = max(results, key=results.get)
    
    print("Likelihood Scores:")
    for diagnosis, score in results.items():
        print(f"{diagnosis}: {score}")
    
    print("\n------------------------------------------------\n")
    print(f"The most likely diagnosis is: {best_diagnosis}")
    print(f"Explanation: {diagnoses[best_diagnosis]['explanation']}")
    print("\nThe scoring equation for the most likely diagnosis is:")
    
    # Print each number in the final equation as requested
    equation_for_best = equations[best_diagnosis]
    final_score_for_best = results[best_diagnosis]
    # Split the string to print each number separately
    numbers_in_equation = equation_for_best.replace(' + ', ' ').replace('- ', '-').split(' ')
    
    # Reconstruct the print output to show each number
    final_equation_str = ""
    for i, num_str in enumerate(numbers_in_equation):
        final_equation_str += num_str
        if i < len(numbers_in_equation) - 1:
            final_equation_str += " + "
    final_equation_str = final_equation_str.replace('+ -', '- ')
    
    print(f"{final_equation_str} = {final_score_for_best}")

diagnose_ankle_pain()