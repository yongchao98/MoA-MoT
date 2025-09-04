import sys
from textwrap import dedent

def check_correctness_of_llm_answer():
    """
    Checks the correctness of the LLM's answer regarding optical isomerism.

    The function analyzes each compound based on chemical principles and, where
    possible, uses the RDKit library to programmatically verify chirality.

    Returns:
        str: "Correct" if the LLM's answer is correct, otherwise a string
             explaining the reason for the error.
    """
    try:
        from rdkit import Chem
        rdkit_available = True
    except ImportError:
        rdkit_available = False
        print("Warning: RDKit not installed. Falling back to logic-based checks only.", file=sys.stderr)

    # The LLM's answer identifies compounds 1 and 2 as optically active, corresponding to option A.
    # We will verify this by checking each compound.
    
    # --- Analysis of each compound ---
    
    # 1. dimethyl 6,6'-dinitro-[1,1'-biphenyl]-2,2'-dicarboxylate
    # This is a classic example of atropisomerism. Chirality arises from restricted
    # rotation around the single bond connecting the two phenyl rings due to bulky
    # ortho substituents (-NO2 and -COOCH3). Since each ring is asymmetrically
    # substituted and the molecule is forced into a non-planar conformation, it
    # lacks a plane of symmetry and is chiral.
    # A simple programmatic check for chiral centers will fail here. We must apply this chemical rule.
    shows_optical_isomerism_1 = True
    
    # 2. methyl 2-hydroxypropanoate (methyl lactate)
    # Structure: CH3-CH(OH)-COOCH3
    # The carbon atom attached to the -OH group is bonded to four different groups:
    # -H, -OH, -CH3, and -COOCH3. It is a chiral center.
    shows_optical_isomerism_2 = True
    if rdkit_available:
        mol2 = Chem.MolFromSmiles("CC(O)C(=O)OC")
        # FindMolChiralCenters returns a list of (atom_idx, chirality_tag)
        # An empty list means no chiral centers were found.
        if not Chem.FindMolChiralCenters(mol2, includeUnassigned=True):
            return dedent("""
                Incorrect.
                Reason: The code's analysis for compound 2 (methyl 2-hydroxypropanoate) contradicts the expected result.
                This molecule has a chiral center at C2 (bonded to H, OH, CH3, and COOCH3) and should be optically active.
                The RDKit check failed to identify this.
            """).strip()

    # 3. benzophenone
    # Structure: C6H5-C(=O)-C6H5
    # This molecule is symmetric. It has a C2 axis and a plane of symmetry.
    # Any molecule with a plane of symmetry is achiral.
    shows_optical_isomerism_3 = False
    if rdkit_available:
        mol3 = Chem.MolFromSmiles("O=C(c1ccccc1)c1ccccc1")
        if Chem.FindMolChiralCenters(mol3, includeUnassigned=True):
            return dedent("""
                Incorrect.
                Reason: The code's analysis for compound 3 (benzophenone) contradicts the expected result.
                This molecule is achiral due to symmetry.
                The RDKit check incorrectly identified a chiral center.
            """).strip()

    # 4. dimethyl fumarate
    # Structure: CH3OOC-CH=CH-COOCH3 (trans isomer)
    # This molecule is planar and has a center of inversion. It is achiral.
    shows_optical_isomerism_4 = False
    if rdkit_available:
        mol4 = Chem.MolFromSmiles("COC(=O)/C=C/C(=O)OC")
        if Chem.FindMolChiralCenters(mol4, includeUnassigned=True):
            return dedent("""
                Incorrect.
                Reason: The code's analysis for compound 4 (dimethyl fumarate) contradicts the expected result.
                This molecule is planar and achiral.
                The RDKit check incorrectly identified a chiral center.
            """).strip()

    # --- Final Verification ---
    # The LLM's answer is A, which means compounds 1 and 2 are the correct ones.
    llm_selection = {1, 2}
    
    our_selection = set()
    if shows_optical_isomerism_1: our_selection.add(1)
    if shows_optical_isomerism_2: our_selection.add(2)
    if shows_optical_isomerism_3: our_selection.add(3)
    if shows_optical_isomerism_4: our_selection.add(4)

    if llm_selection == our_selection:
        return "Correct"
    else:
        return dedent(f"""
            Incorrect.
            Reason: The LLM's answer states that the optically active compounds are {sorted(list(llm_selection))}.
            However, based on chemical principles, the correct set of optically active compounds is {sorted(list(our_selection))}.
            The LLM's analysis of which compounds show optical isomerism is flawed.
        """).strip()

# To run this check, you would need to have rdkit installed:
# pip install rdkit
result = check_correctness_of_llm_answer()
print(result)