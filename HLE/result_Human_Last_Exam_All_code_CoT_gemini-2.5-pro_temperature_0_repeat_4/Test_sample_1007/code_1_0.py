def find_ballet_step():
    """
    Analyzes classical ballet steps to find which one has the same
    starting and ending leg position.
    """
    steps_analysis = {
        'A': {
            'name': 'Entrechat six',
            'description': 'Starts in fifth position, jumps vertically, beats six times, and lands in the same fifth position.'
        },
        'B': {
            'name': 'Échappé battu changé',
            'description': 'Starts in fifth position, jumps to second, then lands back in fifth position but with the feet changed.'
        },
        'C': {
            'name': 'Assemblé',
            'description': 'Starts on one foot, jumps, and lands on two feet in fifth position.'
        },
        'D': {
            'name': 'Glissade derrière',
            'description': 'Starts in fifth position, travels, and closes in fifth position with the feet having changed places.'
        },
        'E': {
            'name': 'Gargouillade',
            'description': 'A complex jump that typically finishes with the feet changed from their starting position.'
        }
    }

    correct_answer = 'A'
    
    print("Analysis of Ballet Steps:")
    for option, details in steps_analysis.items():
        print(f"\nOption {option}: {details['name']}")
        print(f"   - Movement: {details['description']}")

    print("\nConclusion:")
    print("An Entrechat six is a vertical jump that begins and ends in the same position without traveling or changing the feet.")
    print(f"Therefore, the correct option is {correct_answer}.")

find_ballet_step()