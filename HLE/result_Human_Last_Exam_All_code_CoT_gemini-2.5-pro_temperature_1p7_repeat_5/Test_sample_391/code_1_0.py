import sys
import io

# Backup original stdout
original_stdout = sys.stdout
# Redirect stdout to a string stream
sys.stdout = captured_output = io.StringIO()

def solve_antibiotic_choice():
    """
    Analyzes antibiotic susceptibility results to determine the best treatment option.
    """
    culture_results = {
        'Amoxicillin': 'R',
        'Ampicillin/Sulbactam': 'R',
        'Cefalexin': 'R',
        'Cefazolin': 'R',
        'Cefuroxime': 'R',
        'Ciprofloxacin': 'R',
        'Clindamycin': 'S',
        'Erythromycin': 'R',
        'Gentamicin': 'R',
        'Levofloxacin': 'R',
        'Tetracycline': 'I',
        'Vancomycin': 'S',
        'Trimethoprim/Sulfamethoxazole': 'S',
        'Moxifloxacin': 'R',
        'Linezolid': 'S'
    }

    answer_choices = {
        'A': ['Amoxicillin', 'Ciprofloxacin', 'Cefazolin'],
        'B': ['Clindamycin', 'Amoxicillin', 'Tetracycline'],
        'C': ['Vancomycin', 'Linezolid', 'Clindamycin'],
        'D': ['Vancomycin', 'Linezolid', 'Tetracycline'],
        'E': ['Erythromycin', 'Trimethoprim/Sulfamethoxazole', 'Linezolid'],
        'F': ['Gentamicin', 'Levofloxacin', 'Tetracycline'],
        'G': ['Tetracycline', 'Vancomycin', 'Moxifloxacin', 'Trimethoprim/Sulfamethoxazole']
    }

    print("Analyzing treatment options based on culture results:")
    print("S = Susceptible, I = Intermediate, R = Resistant\n")

    best_choice = None
    best_choice_score = -1

    for choice_letter, antibiotics in answer_choices.items():
        is_reasonable = True
        choice_score = 0
        analysis_str = f"Choice {choice_letter}: "
        
        reasons = []
        for antibiotic in antibiotics:
            result = culture_results.get(antibiotic, 'Unknown')
            reasons.append(f"{antibiotic} ({result})")
            if result == 'R':
                is_reasonable = False
            elif result == 'S':
                choice_score += 2 # Prefer S over I
            elif result == 'I':
                choice_score += 1
        
        analysis_str += ", ".join(reasons)
        if is_reasonable:
            analysis_str += " -> This is a reasonable option."
            # The best choice is the one with the highest score (all 'S' preferred)
            if choice_score > best_choice_score:
                best_choice_score = choice_score
                best_choice = choice_letter
        else:
            analysis_str += " -> This is NOT a reasonable option because it includes a resistant antibiotic."
        
        print(analysis_str)

    print("\nConclusion:")
    print(f"The best option is the one containing only Susceptible ('S') antibiotics.")
    print(f"Based on the analysis, Choice {best_choice} is the most appropriate treatment regimen.")


solve_antibiotic_choice()

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()
print(output)
print("<<<C>>>")
