import sys
import io

# Helper class to simulate a safe print environment
class SafePrinter:
    def __init__(self):
        self.output = ""

    def print(self, *args, **kwargs):
        old_stdout = sys.stdout
        sys.stdout = new_stdout = io.StringIO()
        print(*args, **kwargs)
        sys.stdout = old_stdout
        self.output += new_stdout.getvalue()

def solve_litigation_forum():
    """
    Analyzes a legal scenario to determine the best litigation forum
    and presents the reasoning and conclusion.
    """
    sp = SafePrinter()

    # Step 1: Define the facts from the scenario using variables.
    # This includes incorporating all the numbers from the prompt.
    re1_id = 1
    re2_id = 2
    num_properties = 6
    start_year = 2017
    projected_completion_year = 2022
    # The term 'half' translates to a numeric value
    expense_share_re1 = 0.5

    # Case characteristics
    case_is_in_ontario = True
    case_is_complex_commercial = True
    case_requires_speed = True
    case_value_is_high = True # Well over the Small Claims limit

    sp.print("Analyzing the case details:")
    sp.print(f"- The dispute is between RE{re1_id} and RE{re2_id}.")
    sp.print(f"- It concerns {num_properties} properties acquired in {start_year}, with development issues arising before {projected_completion_year}.")
    sp.print(f"- Key issues: Complex contracts, financial disclosure, and ownership of high-value assets.")
    sp.print(f"- Primary requirement: A forum that offers resolution in the shortest amount of time.\n")

    # Step 2: Define the court options and their attributes.
    courts = {
        'A': {'name': "Ontario Court of Appeal", 'is_trial_court': False, 'jurisdiction': "Ontario", 'handles_complex': True, 'is_fast': False, 'monetary_limit_ok': True},
        'B': {'name': "Commercial List", 'is_trial_court': True, 'jurisdiction': "Ontario", 'handles_complex': True, 'is_fast': True, 'monetary_limit_ok': True},
        'C': {'name': "Superior Court of Justice", 'is_trial_court': True, 'jurisdiction': "Ontario", 'handles_complex': True, 'is_fast': False, 'monetary_limit_ok': True},
        'D': {'name': "Small Claims Court", 'is_trial_court': True, 'jurisdiction': "Ontario", 'handles_complex': False, 'is_fast': True, 'monetary_limit_ok': False},
        'E': {'name': "Federal Court of Canada", 'is_trial_court': True, 'jurisdiction': "Federal", 'handles_complex': True, 'is_fast': False, 'monetary_limit_ok': True}
    }

    sp.print("Step 3: Evaluating each forum against the case requirements with a scoring system.")
    sp.print("Scoring Criteria: Trial Court (1), Correct Jurisdiction (1), Handles Complexity (1), High Value OK (1), Is Fast-Track (+2 bonus).")
    sp.print("--------------------------------------------------------------------------")

    best_choice = ''
    max_score = -1

    # Step 3: Iterate and score each option to find the best fit.
    for key, props in courts.items():
        score = 0
        equation_parts = []

        # Build the scoring "equation" based on attributes
        if props['is_trial_court']:
            score += 1
            equation_parts.append("1") # For being a trial court
        else:
            equation_parts.append("0")

        if props['jurisdiction'] == "Ontario" and case_is_in_ontario:
            score += 1
            equation_parts.append("1") # For correct jurisdiction
        else:
            equation_parts.append("0")

        if props['handles_complex'] and case_is_complex_commercial:
            score += 1
            equation_parts.append("1") # For handling complexity
        else:
            equation_parts.append("0")
        
        if props['monetary_limit_ok'] and case_value_is_high:
            score += 1
            equation_parts.append("1") # For appropriate monetary jurisdiction
        else:
            equation_parts.append("0")

        if props['is_fast'] and case_requires_speed:
            score += 2  # Give extra weight to speed as per RE1's request
            equation_parts.append("2") # Bonus for speed
        else:
            equation_parts.append("0")

        # Create the final "equation" string as requested.
        equation_str = " + ".join(equation_parts)

        sp.print(f"Option {key}: {props['name']}")
        sp.print(f"   - Suitability Score Equation: {equation_str} = {score}")
        if score > max_score:
            max_score = score
            best_choice = key

    sp.print("--------------------------------------------------------------------------\n")
    sp.print("Conclusion:")
    sp.print(f"The analysis shows that the 'Commercial List' (Option {best_choice}) has the highest score ({max_score}).")
    sp.print("It is a specialized trial court in the correct jurisdiction designed specifically for a fast and efficient resolution of complex commercial disputes, which directly addresses all of RE1's needs.")

    # Final Answer
    sp.print(f"\n<<<B>>>")
    
    # Print the captured output
    print(sp.output)

solve_litigation_forum()