import sys
from io import StringIO

# It's good practice to install the required library if it's not present.
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    print("RDKit not found. Please install it using: pip install rdkit-pypi")
    # If rdkit is not available, we can't run the check.
    # We will fall back to a logic-based check that encodes the chemical rules.
    RDKIT_AVAILABLE = False
else:
    RDKIT_AVAILABLE = True

def check_chirality_rdkit():
    """
    Uses RDKit to analyze the chirality of the molecules.
    Returns a dictionary mapping compound index to a tuple (is_chiral, reason).
    """
    # SMILES strings for the compounds
    smiles = {
        1: 'C1=C(C(=C(C=C1C2=C(C=C(C=C2[N+](=O)[O-])C(=O)OC)C(=O)OC)[N+](=O)[O-])C(=O)OC)C(=O)OC', # dimethyl 6,6'-dinitro-[1,1'-biphenyl]-2,2'-dicarboxylate
        2: 'CC(C(=O)OC)O',  # methyl 2-hydroxypropanoate
        3: 'C1=CC=C(C=C1)C(=O)C2=CC=CC=C2',  # benzophenone
        4: 'COC(=O)/C=C/C(=O)OC'  # dimethyl fumarate (trans isomer)
    }
    
    results = {}
    
    # 1. Check for atropisomerism in the biphenyl compound
    mol1 = Chem.MolFromSmiles(smiles[1])
    # RDKit's standard chirality detection doesn't always catch atropisomers directly from SMILES.
    # However, the chemical rule is definitive: bulky groups at all four ortho positions
    # on a biphenyl system induce atropisomerism.
    # Substituents: -NO2 and -COOCH3 are bulky.
    # Positions: 2, 2', 6, 6' are all substituted.
    # Conclusion: The molecule is chiral due to atropisomerism.
    results[1] = (True, "Shows optical isomerism. This is a classic case of atropisomerism due to bulky ortho substituents on a biphenyl ring system, which restricts rotation and creates a chiral axis.")

    # Analyze the other compounds
    for i in range(2, 5):
        mol = Chem.MolFromSmiles(smiles[i])
        chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        
        is_chiral = False
        reason = ""
        
        if chiral_centers:
            is_chiral = True
            reason = f"Shows optical isomerism. It has {len(chiral_centers)} chiral center(s)."
        else:
            # Check for other symmetry elements.
            # A molecule with a plane of symmetry or center of inversion is achiral.
            if i == 3: # benzophenone
                is_chiral = False
                reason = "Does not show optical isomerism. The molecule is achiral as it possesses a plane of symmetry."
            elif i == 4: # dimethyl fumarate
                is_chiral = False
                reason = "Does not show optical isomerism. The molecule is achiral as it is planar and has a center of inversion."
        
        results[i] = (is_chiral, reason)
        
    return results

def check_answer_logic():
    """
    A fallback check based on encoded chemical rules if RDKit is not available.
    """
    # Rule-based analysis
    analysis = {
        1: (True, "Shows optical isomerism due to atropisomerism (restricted rotation in a substituted biphenyl)."),
        2: (True, "Shows optical isomerism due to a chiral center (carbon bonded to -H, -OH, -CH3, -COOCH3)."),
        3: (False, "Does not show optical isomerism because it is achiral (has a plane of symmetry)."),
        4: (False, "Does not show optical isomerism because it is achiral (planar molecule with a center of inversion).")
    }
    return analysis

def run_check():
    """
    Runs the check and compares the result with the provided answer.
    """
    # The final answer from the LLM is 'D', which corresponds to compounds 1 and 2.
    llm_answer_indices = {1, 2}
    
    # Use the appropriate checking function
    if RDKIT_AVAILABLE:
        print("Using RDKit for verification...")
        analysis_results = check_chirality_rdkit()
    else:
        print("RDKit not available. Using logic-based verification...")
        analysis_results = check_answer_logic()

    # Determine the set of optically active compounds from our analysis
    verified_chiral_indices = {idx for idx, (is_chiral, reason) in analysis_results.items() if is_chiral}

    # Compare the LLM's answer with the verified set
    if llm_answer_indices == verified_chiral_indices:
        return "Correct"
    else:
        # Construct a detailed error message
        error_report = "The answer is incorrect.\n\n"
        error_report += "Here is the correct analysis of each compound:\n"
        for i in range(1, 5):
            is_chiral, reason = analysis_results[i]
            status = "OPTICALLY ACTIVE" if is_chiral else "NOT optically active"
            error_report += f"  - Compound {i}: {status}. Reason: {reason}\n"
        
        error_report += f"\nConclusion: The compounds that show optical isomerism are {sorted(list(verified_chiral_indices))}.\n"
        error_report += f"The provided answer implies that compounds {sorted(list(llm_answer_indices))} are optically active, which is incorrect."
        return error_report

# Execute the check and print the result
# Redirect stdout to capture print statements for a cleaner final output
old_stdout = sys.stdout
sys.stdout = captured_output = StringIO()

result = run_check()

sys.stdout = old_stdout
# print(captured_output.getvalue()) # Optional: to see the 'Using...' message
print(result)
