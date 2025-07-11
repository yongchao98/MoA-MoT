import textwrap

def analyze_thermodynamics_question():
    """
    Analyzes the multiple-choice question about the limitations of bulk
    calorimetric experiments for nucleic acid thermodynamics.
    """
    
    title = "Analysis of Limitations in Bulk Calorimetric Experiments"
    print(title)
    print("=" * len(title))

    # A dictionary to hold the analysis for each choice
    analysis = {
        'A': "Heat capacity change is assumed to be zero. This is a common simplification in the analytical MODEL, but not a fundamental limitation of the bulk EXPERIMENT itself. The experimental data can be analyzed with more complex models that include a non-zero heat capacity change.",
        'B': "The NNBP parameters are T-independent. This is a consequence of the simplified model mentioned in A. It's a limitation of the data analysis approach, not the experimental measurement, which intrinsically captures effects across a temperature range.",
        'C': "Impossibility to capture heterogeneity in bulk experiments. This is a core, fundamental limitation. Bulk experiments measure the average behavior of an entire population of molecules. If the sample contains different structures or folding pathways (heterogeneity), the experiment averages them all together, masking the properties of individual sub-populations. This is the correct answer.",
        'D': "Temperature oscillations... are too large. This describes an experimental imperfection. The question asks for a limitation that exists even under 'ideal experimental conditions', which would imply perfect temperature control.",
        'E': "Temperature cannot be controlled. This is factually incorrect. Calorimetric melting experiments are defined by precise control and variation of temperature."
    }

    # Print the analysis for each option
    for choice, text in analysis.items():
        # Highlight the correct choice
        prefix = "Correct Choice ->" if choice == 'C' else "Incorrect Choice:"
        wrapped_text = textwrap.fill(text, width=80)
        print(f"\n{prefix} ({choice})")
        print("-" * 20)
        print(wrapped_text)

    print("\n" + "=" * len(title))
    print("Conclusion: The inherent limitation of averaging over a large population of molecules is the inability to see the differences among them.")


# Run the analysis
analyze_thermodynamics_question()