import sys

class ChemicalAnalysis:
    def __init__(self, candidate_id, description, reaction_type):
        self.id = candidate_id
        self.description = description
        self.reaction_type = reaction_type
        self.checks = {}

    def run_checks(self, target):
        # Check 1: Reaction type must be intramolecular to form the fused system directly.
        self.checks['reaction_type'] = self.reaction_type == target['reaction_type']
        if not self.checks['reaction_type']:
            return False

        # Check 2: For IMDA, tether length must be 4 atoms to form a second 6-membered ring.
        self.checks['tether_length'] = self.tether_length == target['tether_length']
        if not self.checks['tether_length']:
            return False

        # Check 3: The connectivity of the final product must match the target.
        # This involves checking the double bond position and substituent arrangement.
        self.checks['product_double_bond'] = self.product_double_bond_position == target['double_bond_position']
        self.checks['product_substituents'] = self.product_substituent_connectivity == target['substituent_connectivity']
        
        return all(self.checks.values())

    def get_failure_reason(self):
        if not self.checks.get('reaction_type', True):
            return f"Candidate {self.id} is an intermolecular reaction, which is not the most direct route to a fused bicyclic system like the target."
        if not self.checks.get('tether_length', True):
            return f"Candidate {self.id} has a tether length of {self.tether_length}, but a 4-atom tether is required to form the bicyclo[4.4.0]decane skeleton."
        if not self.checks.get('product_double_bond', True):
            return f"Candidate {self.id} has a diene at {self.diene_position}, which would form a double bond at {self.product_double_bond_position}. The target requires a double bond at C3=C4."
        if not self.checks.get('product_substituents', True):
            return f"Candidate {self.id} would result in a product with connectivity '{self.product_substituent_connectivity}'. The target requires '{target['substituent_connectivity']}'."
        return "No failure found."

# Define target molecule properties based on its IUPAC name
target = {
    'name': 'methyl 2-propyl-1,2,4a,5,6,7,8,8a-octahydronaphthalene-1-carboxylate',
    'reaction_type': 'intramolecular',
    'tether_length': 4,
    'double_bond_position': 'C3=C4', # In final product numbering
    'substituent_connectivity': 'C1(COOMe)-C2(Propyl)-C3=C4'
}

# Define candidates and their properties for analysis
candidates = []

# Candidate A
cand_A = ChemicalAnalysis('A', 'methyl (2E,4E,10Z)-tetradeca-2,4,10-trienoate', 'intramolecular')
cand_A.diene_position = (2, 3, 4, 5)
cand_A.dienophile_position = (10, 11)
cand_A.tether_length = cand_A.dienophile_position[0] - cand_A.diene_position[3] - 1 # C6, C7, C8, C9 -> 4 atoms
# Product double bond is between central carbons of diene (C3, C4)
cand_A.product_double_bond_position = 'C3=C4'
# Product has ester on C2 of diene, which is part of the new double bond.
cand_A.product_substituent_connectivity = 'C1(COOMe)=C2-...' # Incorrect connectivity
candidates.append(cand_A)

# Candidate B
cand_B = ChemicalAnalysis('B', 'Cyclohexene and methyl 2,3-dimethylenehexanoate', 'intermolecular')
candidates.append(cand_B)

# Candidate C
cand_C = ChemicalAnalysis('C', '1-vinylcyclohex-1-ene and methyl hex-2-ynoate', 'intermolecular')
candidates.append(cand_C)

# Candidate D
cand_D = ChemicalAnalysis('D', 'methyl (2E,8E,10E)-tetradeca-2,8,10-trienoate', 'intramolecular')
cand_D.dienophile_position = (2, 3)
cand_D.diene_position = (8, 9, 10, 11)
cand_D.tether_length = cand_D.diene_position[0] - cand_D.dienophile_position[1] - 1 # C4, C5, C6, C7 -> 4 atoms
# Product double bond is between central carbons of diene (C9, C10).
# Mapping to IUPAC: C1(COOMe) <- C2(prec), C2(Propyl) <- C11(prec), C3 <- C10(prec), C4 <- C9(prec)
# So the double bond is at C3=C4 in the final product.
cand_D.product_double_bond_position = 'C3=C4'
# The connectivity matches the target.
cand_D.product_substituent_connectivity = 'C1(COOMe)-C2(Propyl)-C3=C4'
candidates.append(cand_D)

# --- Main Checking Logic ---
llm_answer = 'D'
correct_candidate_id = None
failure_reasons = {}

for cand in candidates:
    is_correct = cand.run_checks(target)
    if is_correct:
        correct_candidate_id = cand.id
    else:
        failure_reasons[cand.id] = cand.get_failure_reason()

if llm_answer == correct_candidate_id:
    print("Correct")
else:
    if llm_answer in failure_reasons:
        reason = f"The provided answer '{llm_answer}' is incorrect. Reason: {failure_reasons[llm_answer]}"
    else:
        reason = f"The provided answer '{llm_answer}' is incorrect. The correct answer is '{correct_candidate_id}'."
        if correct_candidate_id in failure_reasons:
             reason += f" The provided answer '{llm_answer}' was not identified as incorrect because: {failure_reasons[correct_candidate_id]}"
    print(reason)
