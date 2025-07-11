import sys
import io

# Keep the original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

def solve_agitation_scenario():
    """
    This function models the clinical decision for an agitated patient
    by scoring each treatment option based on safety and effectiveness.
    """
    initial_dose = {'drug': 'olanzapine', 'dose_mg': 5, 'route': 'IM'}

    options = {
        'A': {'description': '2mg IV lorazepam', 'effectiveness': 4, 'safety_route': -5, 'safety_drug': -2},
        'B': {'description': '2mg IM lorazepam + 5mg olanzapine IM', 'effectiveness': 9, 'safety_route': 5, 'safety_drug': -2},
        'C': {'description': 'Verbal de-escalation', 'effectiveness': -10, 'safety_route': 0, 'safety_drug': 0},
        'D': {'description': '10mg IM olanzapine', 'effectiveness': 6, 'safety_route': 5, 'safety_drug': -1},
        'E': {'description': '10mg IM olanzapine + 2mg IM lorazepam', 'effectiveness': 10, 'safety_route': 5, 'safety_drug': -6}
    }

    print("Evaluating treatment options for a patient with acute agitation, unresponsive to 5mg IM olanzapine.")
    print("Scoring based on Effectiveness, Route Safety, and Drug Safety.\n")
    
    best_option = None
    max_score = -float('inf')

    reasoning = {
        'A': "IV route is dangerous in a violent patient and parenteral lorazepam + olanzapine has warnings.",
        'B': "Excellent effectiveness from combination therapy. The total olanzapine dose (5mg initial + 5mg new = 10mg) is appropriate. IM route is safe.",
        'C': "Inappropriate and unsafe as the patient is already violent and pharmacologic treatment has been initiated.",
        'D': "Reasonable to increase the dose of a single drug, but likely less effective than combination therapy for refractory agitation.",
        'E': "Very high risk of over-sedation. The total dose (5mg initial + 10mg new = 15mg olanzapine) plus 2mg lorazepam is excessive for a patient with an unknown history."
    }

    final_equations = {}

    for option, values in options.items():
        score = values['effectiveness'] + values['safety_route'] + values['safety_drug']
        final_equations[option] = f"Option {option} ('{values['description']}'):\nScore Calculation: {values['effectiveness']} (Effectiveness) + {values['safety_route']} (Route Safety) + {values['safety_drug']} (Drug Safety) = {score}\nReasoning: {reasoning[option]}\n"
        if score > max_score:
            max_score = score
            best_option = option

    for option in sorted(final_equations.keys()):
        print(final_equations[option])

    print(f"The best next step is Option {best_option} with a score of {max_score}.")

solve_agitation_scenario()

# Get the content of the buffer
output = captured_output.getvalue()
# Restore the original stdout
sys.stdout = original_stdout

# Print the captured output
print(output)
print("<<<B>>>")