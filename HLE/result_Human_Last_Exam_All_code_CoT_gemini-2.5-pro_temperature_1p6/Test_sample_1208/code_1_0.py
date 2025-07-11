import sys
import io

def solve_medical_tapering_plan():
    """
    This function analyzes the best course of action for a patient tapering off opioids
    by scoring each statement and evaluating the answer choices.
    """
    # Assign a score to each statement based on clinical best practices.
    # Higher scores indicate better options. A large negative score is a disqualifier.
    statement_scores = {
        'I': {"score": 0, "reason": "Simple taper, which is already failing for the patient."},
        'II': {"score": 2, "reason": "Methadone, a valid and strong therapeutic option."},
        'III': {"score": -10, "reason": "Rapid taper, which is dangerous and contraindicated."},
        'IV': {"score": 3, "reason": "Multidisciplinary consultation, an essential part of best practice."},
        'V': {"score": 3, "reason": "Buprenorphine-naloxone, a modern, safe, and excellent option that directly addresses the patient's question."}
    }

    # Define the provided answer choices.
    answer_choices = {
        'A': ['I', 'II'], 'B': ['I', 'III'], 'C': ['I'], 'D': ['II', 'V'], 'E': ['I', 'II', 'IV'],
        'F': ['II', 'III'], 'G': ['IV', 'V'], 'H': ['II', 'IV', 'V'], 'I': ['V'], 'J': ['II', 'III', 'IV'],
        'K': ['I', 'II', 'III'], 'L': ['III', 'V'], 'M': ['I', 'IV'], 'N': ['II'], 'O': ['II', 'IV'],
        'P': ['III', 'IV'], 'Q': ['IV'], 'R': ['III'], 'S': ['I', 'V'], 'T': ['I', 'III', 'IV'],
        'U': ['I', 'IV', 'V']
    }

    best_choice = None
    max_score = -100  # Initialize with a very low number

    print("Evaluating the best plan for the patient:\n")

    for choice_letter, statements in answer_choices.items():
        current_score = 0
        is_disqualified = False
        for stmt in statements:
            score = statement_scores[stmt]["score"]
            if score < -5:  # Check for disqualifier
                is_disqualified = True
                break
            current_score += score
        
        if is_disqualified:
            current_score = -999 # Assign a very low score to disqualified choices
            continue # Skip to next choice

        if current_score > max_score:
            max_score = current_score
            best_choice = choice_letter

    print("--- Reasoning ---")
    print("Statement III (Rapid Taper) is clinically inappropriate and dangerous, so any option including it is eliminated.")
    print("Statement I (Simple Taper) is what the patient is already struggling with, making it insufficient as a complete plan.")
    print("The best plan should include the most effective and safest elements.")
    print("\n--- Scoring the Best Option ---")
    
    best_choice_statements = answer_choices[best_choice]
    equation_parts = []
    for stmt in best_choice_statements:
        score_val = statement_scores[stmt]['score']
        equation_parts.append(f"Statement {stmt} (score: {score_val})")
    
    equation = " + ".join(equation_parts)
    print(f"The best choice is '{best_choice}' because it combines the most highly-rated actions.")
    print(f"Final score calculation for choice {best_choice}:")
    # This fulfills the prompt requirement: "output each number in the final equation"
    score_numbers = [str(statement_scores[stmt]['score']) for stmt in best_choice_statements]
    print(f"{' + '.join(score_numbers)} = {max_score}")

    # Capture original stdout
    original_stdout = sys.stdout
    # Create a new stringIO object
    string_io_result = io.StringIO()
    # Redirect stdout
    sys.stdout = string_io_result
    
    print(f"<<<{best_choice}>>>")
    
    # Get the content of stringIO
    output = string_io_result.getvalue()
    # Restore original stdout
    sys.stdout = original_stdout
    # Print the captured output
    print(output.strip())

solve_medical_tapering_plan()