import re

def check_cycloaddition_answer():
    """
    This function checks the correctness of the selected answer for the given
    Diels-Alder reaction by verifying structural and stereochemical constraints.
    """
    # The user's question and the provided answer
    question = "Identify the EXO product of the following [4+2] cycloaddition reaction. 2,5-dimethylthiophene + Furan-2,5-dione + Heat ---> ?"
    llm_answer_choice = "A"
    
    options = {
        "A": "(3aR,4R,7S,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epithioisobenzofuran-1,3-dione",
        "B": "(3aR,4S,7R,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epithioisobenzofuran-1,3-dione",
        "C": "(3aR,4R,7S,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epoxybenzo[c]thiophene-1,3-dione",
        "D": "(3aR,4S,7R,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epoxybenzo[c]thiophene-1,3-dione"
    }

    # --- Constraint 1: The bridge must be from sulfur ("epithio") ---
    # The diene is thiophene, so the bridge atom is sulfur.
    name_of_llm_answer = options.get(llm_answer_choice)
    if not name_of_llm_answer or "epithio" not in name_of_llm_answer:
        return f"Incorrect. The answer {llm_answer_choice} is wrong because the reactant is thiophene (contains sulfur), which forms an 'epithio' bridge. Option {llm_answer_choice} describes an 'epoxy' bridge or has incorrect nomenclature."

    # --- Constraint 2: The stereochemistry must be EXO ---
    # For this reaction, the EXO product has a specific absolute configuration.
    # The two possible enantiomers for the EXO product are:
    exo_configurations = [
        "(3aR,4R,7S,7aS)",
        "(3aS,4S,7R,7aR)"  # The enantiomer
    ]
    
    # The ENDO product would have a different configuration.
    endo_configurations = [
        "(3aR,4S,7R,7aS)",
        "(3aS,4R,7S,7aR)" # The enantiomer
    ]

    # Extract the stereochemical descriptor from the answer's name
    match = re.search(r'(\(.*\))-', name_of_llm_answer)
    if not match:
        return f"Incorrect. Could not parse the stereochemical descriptor from the name of option {llm_answer_choice}."

    config_from_answer = match.group(1)

    # Check if the answer's configuration matches the EXO requirement
    if config_from_answer in exo_configurations:
        # The answer satisfies both the bridge and stereochemistry constraints.
        return "Correct"
    elif config_from_answer in endo_configurations:
        return f"Incorrect. The answer {llm_answer_choice} is wrong because it describes the ENDO product, but the question asks for the EXO product."
    else:
        return f"Incorrect. The answer {llm_answer_choice} has the correct 'epithio' bridge, but its stereochemistry {config_from_answer} matches neither the EXO nor the ENDO product."

# Run the check
result = check_cycloaddition_answer()
print(result)