import sys

def solve_poem_riddle():
    """
    Analyzes a poem to determine what it describes from a list of choices.
    The analysis is presented as a step-by-step scoring system.
    """
    
    poem_analysis = {
        "clue_1": {
            "text": "Naked, cold",
            "explanation": "This points to a cold, weather-related phenomenon on a bare landscape.",
            "points": { 'A': 1, 'B': 0, 'C': 0, 'D': 0, 'E': 0 }
        },
        "clue_2": {
            "text": "knits a veil... lace and glass",
            "explanation": "This describes creating a delicate, intricate, transparent, and crystalline pattern.",
            "points": { 'A': 1, 'B': 0, 'C': 1, 'D': 0, 'E': 1 }
        },
        "clue_3": {
            "text": "on starwort, grass and meadowsweet",
            "explanation": "This indicates the phenomenon covers plants and flora.",
            "points": { 'A': 1, 'B': 1, 'C': 1, 'D': 0, 'E': 0 }
        },
        "clue_4": {
            "text": "waits for pelted Autumn... to fray",
            "explanation": "The subject ('She') is separate from and will be destroyed by the full force of Autumn (wind/storms).",
            "points": { 'A': 1, 'B': 0, 'C': 0, 'D': -1, 'E': 0 } # Negative point for D as it's a direct contradiction
        }
    }

    options = {
        'A': "The intricate, lace-like patterns of frost during Autumn",
        'B': "A floodplain",
        'C': "A spider spinning her web amongst plants",
        'D': "Autumn as a hunter",
        'E': "A seamstress"
    }

    print("Analyzing the poem by scoring each option against key evidence...\n")

    final_scores = {key: 0 for key in options}
    
    # Store the equation numbers for each option
    equation_parts = {key: [] for key in options}

    for clue_name, clue_data in poem_analysis.items():
        print(f"--- Clue: '{clue_data['text']}' ---")
        print(f"Analysis: {clue_data['explanation']}")
        for option, points in clue_data['points'].items():
            if points != 0:
                print(f"  - Option {option} gets {points} point(s).")
                final_scores[option] += points
                equation_parts[option].append(str(points))
        print("")

    print("--- Final Calculation ---")
    print("Based on the clues, we can form an equation for each option's score.\n")
    
    best_option = ''
    highest_score = -sys.maxsize

    for option_key in sorted(options.keys()):
        score = final_scores[option_key]
        parts = equation_parts[option_key]
        
        # Pad with '0' for clues that didn't award points
        while len(parts) < len(poem_analysis):
            parts.append('0')
            
        equation_str = " + ".join(parts)
        
        print(f"Option {option_key}: {options[option_key]}")
        # The prompt asks to output each number in the final equation.
        print(f"Score Equation: {equation_str} = {score}")
        print("-" * 20)

        if score > highest_score:
            highest_score = score
            best_option = option_key

    print(f"\nConclusion: Option '{best_option}' has the highest score of {highest_score}.")
    print(f"The poem most accurately describes: {options[best_option]}")

solve_poem_riddle()