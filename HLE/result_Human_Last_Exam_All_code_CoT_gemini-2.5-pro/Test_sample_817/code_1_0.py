import textwrap

def analyze_hr4_function():
    """
    This function analyzes the provided multiple-choice options about the plant protein HR4
    and determines the correct answer based on scientific literature.
    """

    # Introduction to the analysis
    print("Analyzing the function of the plant protein HR4...\n")

    # Analysis of each option
    analysis = {
        'A': "It is an interactor of the actin assembly factor ADF3. \n   - Analysis: While the actin cytoskeleton is crucial for plant defense, there is no prominent scientific evidence documenting a direct interaction between HR4 and ADF3. This statement is likely false.",
        'B': "It contributes to the defense against the broad spectrum of powdery mildew pathogens. \n   - Analysis: HR4 does contribute to basal resistance against powdery mildew. However, its effectiveness across a 'broad spectrum' of different powdery mildew species is not its most defining or precisely documented characteristic. This statement is plausible but may be less accurate than a more specific molecular description.",
        'C': "It is targeted to the Extrahaustorial membrane during infection by the powdery mildew pathogen. \n   - Analysis: This is a well-documented finding in plant pathology research. Studies have shown that HR4 protein is specifically recruited to the extrahaustorial membrane (EHM), a unique plant-derived membrane that encases the fungal feeding structure (haustorium). This localization is critical for its role in defense. This statement is strongly supported by evidence.",
        'D': "It regulates the defense modulator PAD4 in plant defense against the Psm. \n   - Analysis: HR4 and PAD4 are both components of the plant immune system, but they function in different, albeit interconnected, pathways. There is no direct evidence that HR4's primary role is to regulate PAD4. This statement is likely false.",
        'E': "HR4 is one of the interactors of the PAD4. \n   - Analysis: Similar to option D, a direct physical interaction between HR4 and PAD4 has not been established as a key mechanism in plant defense. They are not known as primary interactors. This statement is likely false."
    }

    # Print the step-by-step evaluation
    for option, text in analysis.items():
        print(f"--- Evaluating Option {option} ---")
        print(textwrap.fill(text, width=80))
        print("\n")

    # Conclusion
    conclusion = "Based on the analysis, the most accurate and scientifically supported statement is C. The recruitment of HR4 to the extrahaustorial membrane is a key, specific, and well-documented aspect of its function in defending against powdery mildew."
    print("--- Conclusion ---")
    print(textwrap.fill(conclusion, width=80))

    # Final Answer in the required format
    final_answer = "<<<C>>>"
    print(f"\nFinal Answer: {final_answer}")

# Execute the analysis
analyze_hr4_function()