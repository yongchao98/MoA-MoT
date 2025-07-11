import textwrap

def analyze_ftir_data():
    """
    Analyzes FTIR data of tardigrade proteins to determine the structural changes upon gelation.
    """

    # Step 1: Define the standard FTIR peak assignments for protein secondary structures.
    peak_assignments = {
        1652: "Alpha-helix",
        1618: "Beta-sheet",
        1680: "Beta-sheet (anti-parallel)",
        1645: "Disordered / Random Coil"
    }

    print("Step 1: Assigning FTIR peaks to protein secondary structures.")
    for peak, structure in peak_assignments.items():
        print(f"- The peak at {peak} cm^-1 corresponds to {structure}.")
    print("-" * 30)

    # Step 2: Interpret the experimental observations.
    print("Step 2: Interpreting the experimental observations.\n")

    # Concentration Titration
    conc_observation = (
        "Observation from concentration titration: As protein concentration increases, "
        "causing gelation, there is a dual increase in the peaks at 1652 cm^-1 "
        "and 1618 cm^-1."
    )
    conc_interpretation = (
        "Interpretation: The proteins start in a disordered state. Upon gelation, they "
        "simultaneously form both Alpha-helix structures (signal at 1652 cm^-1) and "
        "Beta-sheet structures (signal at 1618 cm^-1)."
    )
    print("--- Concentration Experiment ---")
    print(textwrap.fill(conc_observation, width=80))
    print(textwrap.fill(conc_interpretation, width=80))
    print("\n")

    # Heating Experiment
    heating_observation = (
        "Observation from heating: Upon heating the gel, the peaks at 1618 cm^-1 "
        "and 1680 cm^-1 disappear, while the peak at 1645 cm^-1 grows stronger."
    )
    heating_interpretation = (
        "Interpretation: Heating melts the gel, causing the ordered Beta-sheet "
        "structures (signals at 1618 and 1680 cm^-1) to unfold into a "
        "Disordered / Random Coil state (signal at 1645 cm^-1). This confirms "
        "the gel state is rich in beta-sheets."
    )
    print("--- Heating Experiment ---")
    print(textwrap.fill(heating_observation, width=80))
    print(textwrap.fill(heating_interpretation, width=80))
    print("-" * 30)

    # Step 3: Synthesize and conclude.
    print("Step 3: Synthesizing the results to find the best explanation.\n")
    conclusion = (
        "The evidence from both experiments points to a single conclusion. The proteins "
        "begin as Disordered structures. The process of gelation involves these "
        "structures folding to create a final, ordered state that contains *both* "
        "Alpha-helices and Beta-sheets. Therefore, the correct explanation is the "
        "folding of disordered structures into both beta sheets and alpha helices."
    )
    print(textwrap.fill(conclusion, width=80))

if __name__ == '__main__':
    analyze_ftir_data()
    # Based on the analysis, we select the choice that reflects the folding of
    # disordered structures into both alpha-helices and beta-sheets.
    final_answer = 'I'
    print(f"\n<<<I>>>")