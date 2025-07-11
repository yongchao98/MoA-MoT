import textwrap

def explain_thermodynamics_limitation():
    """
    Analyzes the limitations of bulk nucleic acid melting experiments
    and prints a detailed explanation for the correct answer.
    """
    choices = {
        'A': 'Heat capacity change is assumed to be zero.',
        'B': 'The NNBP parameters are T-independent.',
        'C': 'Impossibility to capture heterogeneity in bulk experiments.',
        'D': 'Temperature oscillations in bulk calorimetry are too large to capture T-dependence.',
        'E': 'Temperature cannot be controlled in calorimetric experiments.'
    }
    
    correct_answer_key = 'C'
    
    print("### Analysis of Limitations in Bulk Nucleic Acid Thermodynamics ###\n")
    
    print(f"The correct choice is: C. {choices[correct_answer_key]}\n")
    
    # Use textwrap for clean printing of long explanation strings
    explanation_wrapper = textwrap.TextWrapper(width=80, initial_indent="  ", subsequent_indent="  ")
    
    print("Detailed Explanation:")
    print("-" * 20)
    
    explanation_c = """A 'bulk' experiment measures the average properties of a vast population of molecules simultaneously. Even with perfect instruments and ideal conditions, this technique has an inherent limitation: it performs an 'ensemble average'. If a sample contains a mixture of different structures (e.g., some fully formed duplexes, some hairpins, some partially melted strands), the bulk measurement (like total heat absorbed) reflects the average behavior of all these states combined. It cannot distinguish or characterize these individual subpopulations. This loss of information about molecular heterogeneity is a fundamental problem that prompted the development of single-molecule techniques."""
    
    print("Why C is correct:")
    print(explanation_wrapper.fill(explanation_c))
    
    print("\nWhy other options are incorrect:")
    
    explanation_ab = """Choices A and B are limitations of the *thermodynamic model* used to interpret the experimental data, not of the *experiment* itself. One can use data from a bulk experiment to fit more complex models where heat capacity change is not zero."""
    print(explanation_wrapper.fill(f"A & B: {explanation_ab}"))
    
    explanation_de = """Choices D and E describe experimental *imperfections* or are factually wrong. The question specifies 'ideal experimental conditions', which rules out instrumental flaws like temperature oscillations. Modern calorimeters offer precise temperature control."""
    print(explanation_wrapper.fill(f"D & E: {explanation_de}"))

# Execute the function to print the analysis
explain_thermodynamics_limitation()