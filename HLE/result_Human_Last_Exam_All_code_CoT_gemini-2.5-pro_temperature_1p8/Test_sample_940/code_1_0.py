import textwrap

def analyze_invasive_species():
    """
    Analyzes a list of species to determine which has had the largest
    negative ecosystem impact as an invasive introduced into New Mexico.
    """
    print("Evaluating which species has the largest negative ecosystem impact in New Mexico...")
    print("-" * 70)

    # Candidate E: Leptinotarsa decemlineata (Colorado potato beetle)
    print("Analysis for Choice E (Leptinotarsa decemlineata):")
    analysis_e = "This species is native to the Southwestern US, including New Mexico. As it is not an 'introduced' species, it is eliminated from consideration based on the prompt's criteria."
    print(textwrap.fill(analysis_e, width=70))
    print("-" * 70)

    # Candidates C and D: Lycorma delicatula and Bombus pascuorum
    print("Analysis for Choices C (Lycorma delicatula) and D (Bombus pascuorum):")
    analysis_cd = "Lycorma delicatula (Spotted Lanternfly) is a major threat in the Eastern US but is not established in New Mexico. Bombus pascuorum (Common Carder Bee) is a European species not considered an invasive species in New Mexico. Since neither is a significant invasive species in the state, they are eliminated."
    print(textwrap.fill(analysis_cd, width=70))
    print("-" * 70)

    # Comparing the remaining candidates
    print("Analysis of remaining candidates A, B, and F:")
    analysis_abf = [
        "A. Apis mellifera (Western Honey Bee): Introduced and vital for agriculture, but can outcompete native pollinators. Its negative impact is often debated and considered less severe than a major disease vector.",
        "F. Maruca vitrata (Bean Pod Borer): An agricultural pest, but its impact is largely confined to specific crops and is not as widespread or severe on the total ecosystem as other candidates.",
        "B. Aedes aegypti (Yellow Fever Mosquito): An introduced invasive species established in southern New Mexico. It is a vector for serious human diseases like Zika, dengue, and chikungunya. The spread of disease has a profound negative impact on human health, a key component of the ecosystem."
    ]
    for item in analysis_abf:
        print(textwrap.fill(item, width=70))
    print("-" * 70)

    # Final Conclusion
    print("Conclusion:")
    conclusion_text = "Comparing the remaining options, the spread of multiple debilitating diseases by Aedes aegypti represents the most severe, direct, and largest negative impact on the ecosystem, including its human inhabitants. Therefore, Aedes aegypti is the correct answer."
    print(textwrap.fill(conclusion_text, width=70))

    final_answer_letter = 'B'
    print(f"\nThe final answer is {final_answer_letter}")

# Run the analysis
analyze_invasive_species()