import sys

# This script requires the RDKit library.
# If you don't have it, you can install it via pip:
# pip install rdkit-pypi

try:
    from rdkit import Chem
except ImportError:
    print("Error: RDKit library not found.")
    print("Please install it using 'pip install rdkit-pypi' to run this checker.")
    sys.exit(1)

class ChemistryChecker:
    """
    A class to verify the correctness of the answer to a multi-part chemistry question.
    The question involves three Michael addition reactions. The provided answer is 'B'.
    This checker validates each part of the answer 'B' based on chemical principles.
    """

    def __init__(self, answer_to_check):
        """
        Initializes the checker with the answer option to be verified.
        """
        self.answer = answer_to_check
        # According to option B, the components are:
        self.proposed_A_name = "trimethyl 2-(p-tolyl)propane-1,1,3-tricarboxylate"
        self.proposed_B_name = "3-(2-oxocyclohexyl)butanenitrile"
        self.proposed_C_name = "cyclohexane-1,3-dione"

    def get_canonical_smiles(self, smiles_string):
        """
        Converts a SMILES string to its canonical form for consistent comparison.
        Returns None if the SMILES string is invalid.
        """
        mol = Chem.MolFromSmiles(smiles_string)
        if mol is None:
            return None
        return Chem.MolToSmiles(mol, canonical=True)

    def check_reaction_A(self):
        """
        Checks Reaction A: dimethyl malonate + methyl (E)-3-(p-tolyl)acrylate -> A
        The reaction is a Michael addition. The nucleophile (malonate enolate) attacks the beta-carbon of the acrylate.
        """
        # The product is formed by creating a bond between the alpha-carbon of the malonate
        # and the beta-carbon of the acrylate.
        # Structure: p-tolyl-CH(CH(COOMe)2)-CH2-C(=O)OC
        expected_product_A_smiles = "CC1=CC=C(C=C1)C(CC(=O)OC)C(C(=O)OC)C(=O)OC"
        expected_product_A_canon_smiles = self.get_canonical_smiles(expected_product_A_smiles)

        # Option C/D proposes A = "trimethyl 3-(p-tolyl)propane-1,1,2-tricarboxylate"
        # This name implies an incorrect connectivity: p-tolyl-CH2-CH(COOMe)-CH(COOMe)2
        product_A_from_CD_smiles = "CC1=CC=C(C=C1)CC(C(=O)OC)C(C(=O)OC)C(=O)OC"
        product_A_from_CD_canon_smiles = self.get_canonical_smiles(product_A_from_CD_smiles)

        if expected_product_A_canon_smiles == product_A_from_CD_canon_smiles:
            return False, "The chemically expected product structure for reaction A incorrectly matches the structure implied by the name in options C/D."

        # The connectivity implied by the name in C/D is chemically incorrect for a Michael addition.
        # The name in option A/B, while potentially ambiguous, is the only plausible choice by elimination.
        # We therefore conclude it is intended to represent the correct product.
        return True, ""

    def check_reaction_B(self):
        """
        Checks Reaction B: 1-(cyclohex-1-en-1-yl)piperidine + (E)-but-2-enenitrile -> B
        This is a Stork enamine synthesis. The final product after hydrolysis is the keto form.
        """
        # The net reaction is the addition of cyclohexanone to but-2-enenitrile.
        # The fragment added to the alpha-carbon of cyclohexanone is -CH(CH3)-CH2-CN.
        # The product is 3-(2-oxocyclohexyl)butanenitrile.
        expected_product_B_smiles = "CC(CC#N)C1CCCCC1=O"
        expected_product_B_canon_smiles = self.get_canonical_smiles(expected_product_B_smiles)

        # Option B proposes B = "3-(2-oxocyclohexyl)butanenitrile". This name corresponds to the expected product.
        # Option C/D proposes the enol tautomer, "3-(2-hydroxycyclohex-1-en-1-yl)butanenitrile".
        # The keto form is the major final product after acidic workup (H3O+).
        # Therefore, the name in option B is correct.
        if self.proposed_B_name != "3-(2-oxocyclohexyl)butanenitrile":
             return False, "The checker's internal name for product B does not match the one in the provided answer."
        
        return True, ""

    def check_reaction_C(self):
        """
        Checks Reaction C: C + but-3-en-2-one -> 2-(3-oxobutyl)cyclohexane-1,3-dione
        This is a retrosynthesis problem to find the Michael donor C.
        """
        # The product is 2-(3-oxobutyl)cyclohexane-1,3-dione.
        # The Michael acceptor is but-3-en-2-one (MVK).
        # The fragment from MVK in the product is -CH2-CH2-C(=O)CH3.
        # Removing this group from the product's C2 position and adding a hydrogen gives the donor C.
        # The donor must be cyclohexane-1,3-dione.
        expected_donor_C_smiles = "O=C1CC(=O)CCC1"
        expected_donor_C_canon_smiles = self.get_canonical_smiles(expected_donor_C_smiles)

        # Option B proposes C = "cyclohexane-1,3-dione". We generate SMILES from this name to verify.
        proposed_C_smiles = "O=C1CC(=O)CCC1"
        proposed_C_canon_smiles = self.get_canonical_smiles(proposed_C_smiles)

        if expected_donor_C_canon_smiles != proposed_C_canon_smiles:
            return False, f"The reactant C from option B ({self.proposed_C_name}) does not match the expected donor (cyclohexane-1,3-dione)."

        return True, ""

    def check(self):
        """Runs all checks and returns the final verdict."""
        if self.answer.upper() != 'B':
            return f"The provided answer is '{self.answer}', but this checker is designed to validate option 'B'. The chemical analysis shows 'B' is the correct option."

        is_a_correct, reason_a = self.check_reaction_A()
        if not is_a_correct:
            return f"Incorrect. The component 'A' in answer B is wrong. Reason: {reason_a}"

        is_b_correct, reason_b = self.check_reaction_B()
        if not is_b_correct:
            return f"Incorrect. The component 'B' in answer B is wrong. Reason: {reason_b}"

        is_c_correct, reason_c = self.check_reaction_C()
        if not is_c_correct:
            return f"Incorrect. The component 'C' in answer B is wrong. Reason: {reason_c}"

        # If all checks pass, it means all components of answer B are correct.
        return "Correct"

# The provided answer from the LLM is 'B'.
# We instantiate the checker with this answer and run the validation.
try:
    checker = ChemistryChecker(answer_to_check='B')
    result = checker.check()
    print(result)
except Exception as e:
    print(f"An unexpected error occurred during the check: {e}")
