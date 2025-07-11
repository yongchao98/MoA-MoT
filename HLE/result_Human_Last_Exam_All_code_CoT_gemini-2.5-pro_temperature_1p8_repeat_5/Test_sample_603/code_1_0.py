def solve_thermodynamics_question():
    """
    Analyzes the limitations of bulk nucleic acid melting experiments
    and identifies the correct answer from a list of choices.
    """

    choices = {
        'A': "Heat capacity change is assumed to be zero.",
        'B': "The NNPB parameters are T-independent.",
        'C': "Impossibility to capture heterogeneity in bulk experiments.",
        'D': "Temperature oscillations in bulk calorimetry are too large to capture T-dependence.",
        'E': "Temperature cannot be controlled in calorimetric experiments."
    }

    correct_answer_key = 'C'

    explanation = """
Explanation:
Bulk experiments measure the average properties of a vast population of molecules. A key limitation of this approach is that it cannot distinguish the behavior of individual molecules or subpopulations within the sample.
This is known as heterogeneity. For nucleic acids, this could include:
- Molecules that have misfolded into alternative, non-native structures.
- The existence of stable intermediate states during the melting process.
- Variations in local environment for different molecules.

A bulk melting curve averages all these different behaviors into a single smooth transition. As a result, the thermodynamic parameters derived are ensemble averages and may not accurately represent any single molecular event. Information about different folding pathways or the existence of distinct subpopulations is lost. Single-molecule techniques (like optical tweezers) were developed specifically to overcome this limitation by observing one molecule at a time.
Therefore, the impossibility of capturing heterogeneity is a fundamental limitation of bulk experiments.
"""

    print(explanation)
    print(f"The correct choice is: {correct_answer_key}")
    
solve_thermodynamics_question()