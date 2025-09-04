def check_drug_discovery_step(selected_option: str) -> str:
    """
    Checks the correctness of the selected option for the bioinformatics question.

    The question asks for the most crucial step before in silico docking for a complex
    molecule with multiple chiral and tautomeric forms.

    Args:
        selected_option: The user's choice, e.g., 'A', 'B', 'C', or 'D'.

    Returns:
        A string indicating "Correct" or explaining why the choice is incorrect.
    """
    correct_answer = 'C'
    reasoning = {
        'A': "Incorrect. While combining in silico and in vitro data is a powerful strategy, it's typically a validation step performed *after* initial computational screening has identified promising candidates. Performing in vitro assays on all possible forms of Xantheraquin *before* any computational prioritization would be impractically expensive and time-consuming.",
        'B': "Incorrect. This is a flawed approach because the most thermodynamically stable form of a molecule in isolation is not necessarily the biologically active form. The specific microenvironment of the protein's binding pocket can stabilize higher-energy tautomers or stereoisomers. Ignoring other forms risks missing the true active compound entirely.",
        'C': "Correct. This is the most crucial and logical preparatory step. Given the complexity of multiple chiral and tautomeric forms, it is essential to first enumerate all relevant possibilities. Then, using computational methods to assess physicochemical properties (like pKa, stability at physiological pH, etc.), one can prioritize the most likely forms for the computationally expensive docking step. This balances comprehensiveness with computational feasibility.",
        'D': "Incorrect. This step is out of sequence. ADME (Absorption, Distribution, Metabolism, Excretion) properties are critical for a drug's overall success, but they are typically evaluated *after* a lead compound with promising target affinity has been identified. The primary goal of the initial docking study is to determine if the molecule binds to the target at all. Focusing on ADME before establishing target affinity is premature."
    }

    selected_option = selected_option.upper()

    if selected_option not in reasoning:
        return f"Invalid option '{selected_option}'. Please choose from A, B, C, or D."

    if selected_option == correct_answer:
        return "Correct"
    else:
        return reasoning[selected_option]

# Example usage:
# print(check_drug_discovery_step('A'))
# print(check_drug_discovery_step('B'))
# print(check_drug_discovery_step('C'))
# print(check_drug_discovery_step('D'))
