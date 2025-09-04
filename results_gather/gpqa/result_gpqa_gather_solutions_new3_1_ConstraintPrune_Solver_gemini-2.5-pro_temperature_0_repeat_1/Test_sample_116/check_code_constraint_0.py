import sys
from io import StringIO

# A check to ensure the rdkit library is installed, which is necessary for this chemical analysis.
try:
    from rdkit import Chem
    from rdkit.Chem.rdmolops import GetPotentialStereo, StereoType
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

def check_answer_correctness():
    """
    Checks the correctness of the answer regarding optical isomerism using the rdkit library.
    
    The function analyzes each molecule for sources of chirality:
    1.  Chiral Centers: Atoms (like carbon) bonded to four different substituents.
    2.  Atropisomerism: Axial chirality arising from restricted rotation around a single bond,
        common in bulky ortho-substituted biphenyls.
    
    It then compares the analysis results with the provided final answer.
    """
    if not RDKIT_AVAILABLE:
        return "Could not perform check: The 'rdkit' library is not installed. Please install it (e.g., 'pip install rdkit') to run this verification code."

    # The final answer from the LLM to be checked.
    llm_answer = "D"

    # Define the compounds using their SMILES representation, a standard way to encode molecular structures.
    compounds = {
        1: ("dimethyl 6,6'-dinitro-[1,1'-biphenyl]-2,2'-dicarboxylate", "COC(=O)c1c(cccc1[N+](=O)[O-])c2c(C(=O)OC)cccc2[N+](=O)[O-]"),
        2: ("methyl 2-hydroxypropanoate", "CC(C(=O)OC)O"),
        3: ("benzophenone", "O=C(c1ccccc1)c2ccccc2"),
        4: ("dimethyl fumarate", "COC(=O)/C=C/C(=O)OC")
    }

    # Map the multiple-choice options to the compounds they represent.
    option_map = {
        "A": {2, 3},
        "B": {3, 4},
        "C": {1, 2, 4},
        "D": {1, 2}
    }

    # Perform chemical analysis for each compound.
    analysis_results = {}
    for num, (name, smiles) in compounds.items():
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            # This case should not happen with valid SMILES.
            analysis_results[num] = (False, f"Error: Could not parse SMILES for {name}")
            continue

        is_chiral = False
        reason = "Achiral (possesses symmetry, no chiral elements found)."

        # 1. Check for atom-centered chirality (e.g., a carbon with 4 different groups).
        #    'includeUnassigned=True' is crucial for finding potential chiral centers.
        chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        if len(chiral_centers) > 0:
            is_chiral = True
            reason = f"Has {len(chiral_centers)} chiral center(s)."
        
        # 2. Check for axial chirality (atropisomerism), which is relevant for compound 1.
        if not is_chiral:
            # GetPotentialStereo identifies potential stereogenic elements, including atropisomers.
            potential_stereo = GetPotentialStereo(mol)
            for stereo_info in potential_stereo:
                if stereo_info.type == StereoType.Bond_Atropisomer:
                    is_chiral = True
                    reason = "Exhibits atropisomerism (a form of axial chirality)."
                    break
        
        analysis_results[num] = (is_chiral, reason)

    # Determine the set of optically active compounds based on our analysis.
    active_compounds_from_analysis = {num for num, (is_active, _) in analysis_results.items() if is_active}

    # Get the set of compounds claimed by the LLM's answer.
    if llm_answer not in option_map:
        return f"Invalid answer format: '{llm_answer}'. Expected one of {list(option_map.keys())}."
    
    active_compounds_from_llm_answer = option_map[llm_answer]

    # Compare the analysis with the LLM's answer.
    if active_compounds_from_analysis == active_compounds_from_llm_answer:
        return "Correct"
    else:
        # If incorrect, provide a detailed explanation.
        error_message = (
            f"Incorrect. The provided answer '{llm_answer}' claims that compounds {sorted(list(active_compounds_from_llm_answer))} are optically active, "
            f"but a chemical analysis shows that compounds {sorted(list(active_compounds_from_analysis))} are the correct ones.\n\n"
            "Here is the detailed breakdown:\n"
        )
        for num, (is_active, reason) in analysis_results.items():
            name = compounds[num][0]
            status = "OPTICALLY ACTIVE" if is_active else "NOT optically active"
            error_message += f"- Compound {num} ({name}): {status}. Reason: {reason}\n"
        
        return error_message.strip()

# Run the check and print the result.
print(check_answer_correctness())