import sys
import io

# A helper class to redirect stdout to capture print statements
class Capturing(list):
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = io.StringIO()
        return self
    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        del self._stringio    # free up some memory
        sys.stdout = self._stdout

def solve_clinical_scenario():
    """
    This function analyzes a clinical scenario by scoring different treatment options
    to determine the best course of action.
    """
    print("Analyzing the clinical statements to find the best approach for the patient...")

    # Assign scores to each statement based on clinical appropriateness.
    # Higher scores indicate better options. Negative scores indicate harmful options.
    statement_scores = {
        'I': 1,   # Maintain current regimen: Insufficient as patient is already facing challenges.
        'II': 5,  # Transition to methadone: A valid and effective option.
        'III': -10,# Initiate rapid taper: Dangerous and contraindicated.
        'IV': 10, # Arrange multidisciplinary consultation: Essential for complex cases, best practice.
        'V': 8    # Prescribe buprenorphine-naloxone: Excellent, safe option, addresses patient's question.
    }
    
    # Define the answer choices as combinations of the statements.
    answer_choices = {
        'A': ['I', 'II'], 'B': ['I', 'III'], 'C': ['I'], 'D': ['II', 'V'],
        'E': ['I', 'II', 'IV'], 'F': ['II', 'III'], 'G': ['IV', 'V'], 'H': ['II', 'IV', 'V'],
        'I': ['V'], 'J': ['II', 'III', 'IV'], 'K': ['I', 'II', 'III'], 'L': ['III', 'V'],
        'M': ['I', 'IV'], 'N': ['II'], 'O': ['II', 'IV'], 'P': ['III', 'IV'],
        'Q': ['IV'], 'R': ['III'], 'S': ['I', 'V'], 'T': ['I', 'III', 'IV'],
        'U': ['I', 'IV', 'V']
    }

    best_choice = None
    max_score = -float('inf')
    
    # Calculate the score for each answer choice
    for choice, statements in answer_choices.items():
        current_score = sum(statement_scores[s] for s in statements)
        if current_score > max_score:
            max_score = current_score
            best_choice = choice

    # Output the result and the calculation for the best choice
    print(f"\nBest choice identified: {best_choice} with a total score of {max_score}.")
    
    print("\nCalculation for the best choice:")
    winning_statements = answer_choices[best_choice]
    equation_parts = [f"Score('{s}')" for s in winning_statements]
    score_values = [str(statement_scores[s]) for s in winning_statements]
    
    print(f"Score('{best_choice}') = {' + '.join(equation_parts)}")
    print(f"             = {' + '.join(score_values)} = {max_score}")

    # The rationale is that the best plan is comprehensive: it involves the best practice process (IV: multidisciplinary team)
    # and considers the top evidence-based treatments (II: methadone, V: buprenorphine). This combination
    # yields the highest score.

    final_answer = f"<<<{best_choice}>>>"
    
    # Capture the output and print it, then print the final answer format separately
    return final_answer

# Run the analysis
with Capturing() as output:
    final_answer = solve_clinical_scenario()

# Print the captured output from the function
for line in output:
    print(line)

# Print the final answer in the required format
print(final_answer)
