import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

def solve_antioxidant_question():
    """
    Identifies and prints the correct answer to the biological question
    based on scientific findings.
    """
    answer_choices = {
        'A': 'Liposoluble antioxidants',
        'B': 'Hydrosoluble antioxidants',
        'C': 'Enzymatic antioxidants',
        'D': 'Photosynthetic pigments',
        'E': 'UV-protective compounds'
    }

    # Based on Shao et al. (2018), the primary initial response involves
    # hydrosoluble antioxidants for this specific strain and condition.
    correct_answer_key = 'B'

    print("The initially activated antioxidants in Microcystis aeruginosa CAAT 2005-3 in response to high temperature exposure (29ºC) are:")
    print(f"Answer: {correct_answer_key}. {answer_choices[correct_answer_key]}")

solve_antioxidant_question()

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()
# Print the captured output to the actual console
print(output)