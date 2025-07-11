import textwrap

def explain_thermodynamics_limitation():
    """
    Analyzes the choices to find the fundamental limitation of bulk melting experiments
    for nucleic acid thermodynamics.
    """
    question = "What is a limitation of bulk calorimetric melting experiments for NA thermodynamics, even under ideal experimental conditions?"

    choices = {
        'A': "Heat capacity change is assumed to be zero.",
        'B': "The NNPB parameters are T-independent.",
        'C': "Impossibility to capture heterogeneity in bulk experiments.",
        'D': "Temperature oscillations in bulk calorimetry are too large to capture T-dependence.",
        'E': "Temperature cannot be controlled in calorimetric experiments."
    }

    print("Analyzing the question: " + question + "\n")
    print("-" * 70)

    # Step 1: Analyze choices related to the thermodynamic model (A & B)
    print("Step 1: Analyzing choices A and B.\n")
    analysis_ab = """These choices refer to assumptions often made when analyzing the data (e.g., in a simple van't Hoff plot), not a fundamental limitation of the bulk measurement technique itself. More complex models that account for heat capacity changes (and thus T-dependent parameters) can be applied to bulk data. So, these are limitations of a simplified analysis, not the experiment itself."""
    print(textwrap.fill(analysis_ab, width=70))
    print("\nVerdict: A and B are incorrect.\n")
    print("-" * 70)

    # Step 2: Analyze choices related to experimental artifacts (D & E)
    print("Step 2: Analyzing choices D and E.\n")
    analysis_de = """Choice E is factually incorrect; the precise control of temperature is the foundation of these experiments. Choice D is also incorrect, as modern calorimeters have very high-precision temperature control. These are not fundamental limitations."""
    print(textwrap.fill(analysis_de, width=70))
    print("\nVerdict: D and E are incorrect.\n")
    print("-" * 70)

    # Step 3: Analyze the choice related to the nature of bulk measurements (C)
    print("Step 3: Analyzing choice C.\n")
    analysis_c = """This choice addresses the core concept of a 'bulk' or 'ensemble' measurement. The experiment measures the average properties of a massive population of molecules. At any point during melting, this population is heterogeneous—it contains a mix of fully folded, partially unfolded, and fully unfolded molecules. The bulk experiment averages over all these distinct molecular states, producing a single smooth curve. It fundamentally cannot resolve the behavior of individual molecules or distinguish between different melting pathways and intermediates. This loss of information about the underlying molecular heterogeneity is a primary limitation that motivated the development of single-molecule techniques."""
    print(textwrap.fill(analysis_c, width=70))
    print("\nVerdict: C is correct.\n")
    print("-" * 70)

    # Final Conclusion
    print("Conclusion: The most fundamental limitation is the averaging over a heterogeneous population of molecules.")
    print(f"The correct answer is C: {choices['C']}")

if __name__ == '__main__':
    explain_thermodynamics_limitation()