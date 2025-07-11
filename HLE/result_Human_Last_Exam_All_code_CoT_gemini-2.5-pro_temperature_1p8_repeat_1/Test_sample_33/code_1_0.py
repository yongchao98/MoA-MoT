import collections

def analyze_ftir_data():
    """
    Analyzes FTIR data of a tardigrade protein to determine structural changes upon gelation.
    """
    
    # Step 1: Define known FTIR peak assignments for protein secondary structures (Amide I band).
    structure_peaks = {
        "Alpha-helix": "around 1650-1658 cm^-1",
        "Beta-sheet (anti-parallel)": "strong peak at 1610-1640 cm^-1 and a weaker peak at 1680-1690 cm^-1",
        "Disordered / Random Coil": "around 1640-1645 cm^-1"
    }

    # Step 2: State the key experimental observations from the concentration titration leading to gelation.
    # The protein starts as disordered and becomes ordered as concentration increases.
    gelation_observations = {
        "peak_increase_1": 1652,  # in cm^-1
        "peak_increase_2": 1618   # in cm^-1
    }
    
    # Initial state of the protein is disordered.
    initial_state = "Disordered structure"
    
    # Step 3: Interpret the observations based on the peak assignments.
    conclusion_parts = []
    
    # Check for alpha-helix formation
    observed_peak_1 = gelation_observations["peak_increase_1"]
    if 1650 <= observed_peak_1 <= 1658:
        conclusion_parts.append(f"the increase in the {observed_peak_1} cm^-1 peak indicates the formation of alpha-helices")

    # Check for beta-sheet formation
    observed_peak_2 = gelation_observations["peak_increase_2"]
    if 1610 <= observed_peak_2 <= 1640:
        conclusion_parts.append(f"the increase in the {observed_peak_2} cm^-1 peak indicates the formation of beta-sheets")
        # The heating experiment confirms this by showing the disappearance of the characteristic anti-parallel
        # beta-sheet peaks at 1618 cm^-1 and 1680 cm^-1 upon unfolding.

    # Step 4: Synthesize the final explanation.
    print("Explanation of the protein behavior:\n")
    print(f"The protein is initially in a '{initial_state}'.")
    print("Upon increasing concentration to form a hydrogel:")
    
    explanation = f"The experimental data shows that {conclusion_parts[0]} and {conclusion_parts[1]}."
    
    print(explanation)
    print("\nTherefore, the overall process is the folding of the initially disordered protein into a structure containing both alpha-helices and beta-sheets.")
    print("\nThis matches the following choice:")
    print("I. Disordered structures fold into beta sheets and alpha helices upon gelation")


# Run the analysis
analyze_ftir_data()

# Final Answer determined by the analysis
print("\n<<<I>>>")