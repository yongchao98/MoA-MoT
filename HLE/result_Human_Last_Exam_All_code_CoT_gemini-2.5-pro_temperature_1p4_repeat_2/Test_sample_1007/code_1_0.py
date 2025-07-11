def analyze_ballet_steps():
    """
    Analyzes classical ballet steps to determine which one ends in the same
    leg position as its starting position.
    """
    steps = {
        "A. Entrechat six": {
            "starts_same_as_ends": False,
            "description": "A jump with six crossings. Starting in 5th with the right foot front, the dancer lands in 5th with the left foot front."
        },
        "B. Échappé battu changé": {
            "starts_same_as_ends": False,
            "description": "A jump from 5th to 2nd position, then a jump back to 5th, changing which foot is in front. The name 'changé' means 'changed'."
        },
        "C. Assemblé": {
            "starts_same_as_ends": True,
            "description": "A jump from one foot to two. The working leg brushes out and is then 'assembled' with the supporting leg to land in the original 5th position."
        },
        "D. Glissade derrière": {
            "starts_same_as_ends": False,
            "description": "A gliding step where the starting front foot closes in the back, changing the front foot."
        },
        "E. Gargouillade": {
            "starts_same_as_ends": False,
            "description": "A complex jump that results in landing with the opposite foot in front of where it started."
        }
    }

    correct_answer_choice = None
    print("Analysis of Ballet Steps:")
    print("=" * 30)

    for choice, details in steps.items():
        print(f"Step: {choice}")
        print(f"Description: {details['description']}")
        if details['starts_same_as_ends']:
            print("Result: The ending and starting positions are the SAME.\n")
            correct_answer_choice = choice
        else:
            print("Result: The ending and starting positions are DIFFERENT.\n")

    print("=" * 30)
    print(f"Conclusion: The correct step is '{correct_answer_choice}', as it is the only one that, in its basic form, finishes in the same position from which it started.")

if __name__ == '__main__':
    analyze_ballet_steps()