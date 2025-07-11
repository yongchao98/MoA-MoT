import sys

def analyze_ftir_data():
    """
    Analyzes FTIR data of tardigrade protein gelation to determine the structural transition.
    """

    # Step 1: Define the relationship between FTIR Amide I peaks and protein secondary structures.
    structure_map = {
        1645: "Disordered / Random Coil (broad peak)",
        1652: "Alpha-Helix",
        1618: "Beta-Sheet",
        1680: "Beta-Sheet (specifically, a component of anti-parallel sheets)"
    }

    print("Step 1: Correlating FTIR peak values to protein secondary structures.")
    print("-" * 60)
    for peak, structure in structure_map.items():
        print(f"A peak around {peak} cm⁻¹ corresponds to: {structure}")
    print("-" * 60)
    print("\n")

    # Step 2: Analyze the observations from the concentration titration experiment.
    print("Step 2: Analyzing the gelation process (concentration increase).")
    print("-" * 60)
    peak1 = 1652
    peak2 = 1618
    print(f"Observation: As concentration increases, there is a dual increase in the peaks at {peak1} cm⁻¹ and {peak2} cm⁻¹.")
    print(f"This means that during gelation, the amounts of both '{structure_map[peak1]}' and '{structure_map[peak2]}' are increasing.")
    print("This strongly suggests that the initial disordered proteins are folding into a mix of both alpha-helices and beta-sheets.")
    print("-" * 60)
    print("\n")
    
    # Step 3: Use the heating experiment to confirm the assignments.
    print("Step 3: Confirming with the heating experiment (gel melting).")
    print("-" * 60)
    disordered_peak = 1645
    sheet_peak1 = 1618
    sheet_peak2 = 1680
    print(f"Observation: Upon heating, the '{structure_map[disordered_peak]}' peak at {disordered_peak} cm⁻¹ grows,")
    print(f"while the '{structure_map[sheet_peak1]}' and '{structure_map[sheet_peak2]}' peaks at {sheet_peak1} cm⁻¹ and {sheet_peak2} cm⁻¹ disappear.")
    print("This confirms that the ordered structures (beta-sheets) unfold back into a disordered state upon heating, which is expected behavior for gel melting.")
    print("-" * 60)
    print("\n")

    # Step 4: Formulate the final conclusion.
    print("Final Conclusion:")
    print("-" * 60)
    print("The evidence shows a transition from an initial disordered state to a final gelled state containing both alpha-helical and beta-sheet structures.")
    print("Therefore, the correct explanation is: Disordered structures fold into beta sheets and alpha helices upon gelation.")
    print("-" * 60)


if __name__ == '__main__':
    analyze_ftir_data()
