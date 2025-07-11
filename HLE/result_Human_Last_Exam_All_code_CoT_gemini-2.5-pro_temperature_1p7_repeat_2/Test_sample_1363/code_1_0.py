def solve_dance_question():
    """
    This script analyzes a technical ballroom dance question and provides the correct answer with reasoning.
    """
    question = "In which dance is it impossible to overturn a reverse turn without disregarding the technique?"
    
    options = {
        'A': {'name': 'Viennese Waltz', 'number': 1},
        'B': {'name': 'English Waltz', 'number': 2},
        'C': {'name': 'European Tango', 'number': 3},
        'D': {'name': 'Slow Foxtrott', 'number': 4},
        'E': {'name': 'Quickstep', 'number': 5}
    }

    correct_choice_letter = 'C'
    correct_choice_info = options[correct_choice_letter]

    # Reasoning:
    # The European Tango is fundamentally different from the other four standard dances.
    # It is a staccato, non-swing dance with no rise & fall or body sway.
    # Turns in Tango are sharp and precise. The Reverse Turn in Tango has a very
    # specific structure. Attempting to "overturn" it (add more rotation) would
    # require introducing continuous movement and sway, which is a direct violation of
    # the core Tango technique. In contrast, the other four are swing dances where
    # overturning figures is a common and accepted practice in advanced choreography.
    
    print("Let the multiple-choice options be represented by numbers:")
    for key, value in options.items():
        print(f"{key} = {value['number']} ({value['name']})")
    
    print("\nThe solution can be found in the following equation, where the number represents the correct choice:")
    
    # The final equation as requested, printing each number
    print(f"Correct Answer = {correct_choice_info['number']}")

    print(f"\nThus, the answer is {correct_choice_letter}: {correct_choice_info['name']}.")

solve_dance_question()