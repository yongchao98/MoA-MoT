import textwrap

def explain_protein_folding():
    """
    Analyzes FTIR data for tardigrade hydrogel proteins and explains the structural transitions.
    """
    # Define the FTIR peaks and their corresponding secondary structures
    peak_assignments = {
        1645: "Random Coil / Disordered Structure (indicated by broad peak)",
        1652: "Alpha Helix",
        1618: "Intermolecular Beta Sheet (often related to aggregation)",
        1680: "Anti-parallel Beta Sheet (often paired with a low-frequency peak like 1618 cm^-1)"
    }

    # Use textwrap for clean, indented output
    wrapper = textwrap.TextWrapper(initial_indent="  - ", subsequent_indent="    ")

    print("Step-by-Step Analysis of Protein Folding Behavior:")
    print("-" * 50)

    # 1. Analyze the initial state
    print("1. Initial State Analysis:")
    print("  The protein is initially mostly disordered, as shown by the strong, broad peak at 1645 cm^-1.")
    print("  The presence of shoulders at 1652, 1618, and 1680 cm^-1 suggests small, pre-existing populations of alpha-helices and beta-sheets.")
    print("\n  Peak Assignments:")
    for peak, assignment in peak_assignments.items():
        print(f"  - {peak} cm^-1: {assignment}")
    
    print("\n" + "-" * 50)

    # 2. Analyze the heating experiment
    print("2. Heating Experiment (Unfolding):")
    text_heating = (f"Upon heating, the peak at 1645 cm^-1 (Disordered) grows stronger, while the peaks at "
                    f"1618 cm^-1 and 1680 cm^-1 (Beta Sheet) disappear. This indicates that the ordered beta-sheet "
                    f"structures are unfolding into random coils due to thermal denaturation. This confirms our peak assignments.")
    print(textwrap.fill(text_heating, width=80, initial_indent="  ", subsequent_indent="  "))

    print("\n" + "-" * 50)

    # 3. Analyze the concentration experiment (Gelation/Folding)
    print("3. Concentration Titration (Gelation):")
    text_concentration = (f"As concentration increases, the hydrogel forms. The data shows a dual increase in the peaks at "
                          f"1652 cm^-1 (Alpha Helix) and 1618 cm^-1 (Beta Sheet). This is the key observation.")
    print(textwrap.fill(text_concentration, width=80, initial_indent="  ", subsequent_indent="  "))
    
    print("\n" + "-" * 50)
    
    # 4. Conclusion
    print("4. Conclusion:")
    conclusion = ("The most logical explanation is that the initially disordered proteins fold and assemble into a "
                  "more complex, ordered structure to form the hydrogel. Since signals for BOTH alpha-helices (1652 cm^-1) "
                  "and beta-sheets (1618 cm^-1) increase during this process, the proteins are folding into a "
                  "conformation containing both types of secondary structures.")
    print(textwrap.fill(conclusion, width=80, initial_indent="  ", subsequent_indent="  "))
    print("\n  This corresponds to answer choice I.")


if __name__ == '__main__':
    explain_protein_folding()
    # The final answer is determined by the logic presented.
    # The analysis shows disordered proteins folding into both alpha-helices and beta-sheets.
    final_answer = 'I'
    print(f"\nFinal Answer: <<<__{final_answer}__>>>")
