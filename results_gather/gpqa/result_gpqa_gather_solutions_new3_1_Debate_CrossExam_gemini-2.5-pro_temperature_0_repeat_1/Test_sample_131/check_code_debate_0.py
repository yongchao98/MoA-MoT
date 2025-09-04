import collections
from math import gcd

def check_nmr_answer():
    """
    This function verifies the correct answer for the NMR spectroscopy problem.
    It models the expected NMR spectrum for each compound, simulates the spectrum
    for a 1:1 mixture for each option, and compares it against the experimental data.
    """

    # Step 1: Define the expected 1H NMR data for each compound based on symmetry.
    # Signal format: (integration, multiplicity)
    compound_data = {
        '1,2,4,5-tetramethylbenzene': {
            'aromatic': [(2, 'singlet')],
            'alkyl': [(12, 'singlet')]
        },
        '1,2,3,5-tetramethylbenzene': {
            # No symmetry -> 2 unique aromatic H, 4 unique methyl groups
            'aromatic': [(1, 'singlet'), (1, 'singlet')],
            'alkyl': [(3, 'singlet'), (3, 'singlet'), (3, 'singlet'), (3, 'singlet')]
        },
        '1,2,3,4-tetramethylbenzene': {
            # Plane of symmetry -> 1 set of aromatic H, 2 sets of methyl groups
            'aromatic': [(2, 'singlet')],
            'alkyl': [(6, 'singlet'), (6, 'singlet')]
        },
        '1,4-diethylbenzene': {
            # High symmetry -> 1 set of aromatic H, 1 set of ethyl groups
            'aromatic': [(4, 'singlet')],
            'alkyl': [(4, 'quartet'), (6, 'triplet')]
        }
    }

    # Step 2: Define the options from the question.
    options = {
        'A': ['1,2,4,5-tetramethylbenzene', '1,2,3,4-tetramethylbenzene'],
        'B': ['1,2,3,5-tetramethylbenzene', '1,4-diethylbenzene'],
        'C': ['1,2,3,4-tetramethylbenzene', '1,2,3,5-tetramethylbenzene'],
        'D': ['1,2,4,5-tetramethylbenzene', '1,2,3,5-tetramethylbenzene']
    }

    # Step 3: Define the experimental data from the question.
    experimental_data = {
        'aromatic': {'count': 2, 'multiplicity': 'singlet', 'ratio': [1, 1]},
        'alkyl': {'count': 3, 'multiplicity': 'singlet', 'ratio': [2, 1, 1]}
    }
    
    # The final answer provided by the LLM to be checked.
    llm_answer = 'A'

    # Helper function to simplify a list of numbers into its simplest integer ratio.
    def simplify_ratio(numbers):
        if not numbers:
            return []
        # Find the greatest common divisor (GCD) of all numbers.
        common_divisor = numbers[0]
        for num in numbers[1:]:
            common_divisor = gcd(common_divisor, num)
        # Divide each number by the GCD and sort for consistent comparison.
        return sorted([num // common_divisor for num in numbers], reverse=True)

    # Step 4: Analyze each option against the experimental data.
    correct_option = None
    analysis_log = {}

    for option_key, compounds in options.items():
        comp1_name, comp2_name = compounds
        comp1_data = compound_data[comp1_name]
        comp2_data = compound_data[comp2_name]

        # Combine signals for a 1:1 mixture.
        mixture_aromatic = comp1_data['aromatic'] + comp2_data['aromatic']
        mixture_alkyl = comp1_data['alkyl'] + comp2_data['alkyl']

        # --- Check Constraints ---
        
        # Constraint 1: All alkyl signals must be singlets.
        if not all(sig[1] == 'singlet' for sig in mixture_alkyl):
            analysis_log[option_key] = "Fails: Contains non-singlet alkyl signals."
            continue

        # Constraint 2: The number of signals must match.
        # This check assumes no accidental overlap, which is valid for the correct answer.
        if len(mixture_aromatic) != experimental_data['aromatic']['count']:
            analysis_log[option_key] = f"Fails: Expected {experimental_data['aromatic']['count']} aromatic signals, but mixture has {len(mixture_aromatic)}."
            continue
        if len(mixture_alkyl) != experimental_data['alkyl']['count']:
            analysis_log[option_key] = f"Fails: Expected {experimental_data['alkyl']['count']} alkyl signals, but mixture has {len(mixture_alkyl)}."
            continue

        # Constraint 3: The integration ratios must match.
        alkyl_integrations = [sig[0] for sig in mixture_alkyl]
        calculated_alkyl_ratio = simplify_ratio(alkyl_integrations)
        expected_alkyl_ratio = sorted(experimental_data['alkyl']['ratio'], reverse=True)
        if calculated_alkyl_ratio != expected_alkyl_ratio:
            analysis_log[option_key] = f"Fails: Expected alkyl ratio {expected_alkyl_ratio}, but got {calculated_alkyl_ratio}."
            continue

        aromatic_integrations = [sig[0] for sig in mixture_aromatic]
        calculated_aromatic_ratio = simplify_ratio(aromatic_integrations)
        expected_aromatic_ratio = sorted(experimental_data['aromatic']['ratio'], reverse=True)
        if calculated_aromatic_ratio != expected_aromatic_ratio:
            analysis_log[option_key] = f"Fails: Expected aromatic ratio {expected_aromatic_ratio}, but got {calculated_aromatic_ratio}."
            continue

        # If all checks pass, this is the correct option.
        correct_option = option_key
        analysis_log[option_key] = "Passes all constraints."

    # Step 5: Compare the code's conclusion with the LLM's answer.
    if correct_option == llm_answer:
        return "Correct"
    elif correct_option is None:
        return f"Incorrect. The provided answer is {llm_answer}, but my analysis shows that no option perfectly matches the criteria. The analysis for each option is: {analysis_log}"
    else:
        return f"Incorrect. The provided answer is {llm_answer}, but the correct answer should be {correct_option}. The analysis for option {llm_answer} shows: {analysis_log.get(llm_answer, 'No analysis found.')}. The correct option {correct_option} passes because: {analysis_log.get(correct_option)}."

# Execute the check and print the result.
result = check_nmr_answer()
print(result)