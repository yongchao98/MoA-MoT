import re
from rdkit import Chem

def check_correctness_of_llm_answer():
    """
    Checks the correctness of the LLM's answer and reasoning for the multi-step synthesis problem.

    The function verifies the following:
    1.  The fundamental chemical constraint that the final product must have cis-dicarboxylate groups,
        originating from maleic anhydride.
    2.  The LLM's specific claims about which isomers are cis/trans.
    3.  The LLM's claims about the identity of the options.
    4.  The correctness of the LLM's embedded python code output.
    """
    options = {
        'A': {
            "name": "dimethyl (1S,4R,4aR,4bR,5S,8R,8aS,8bS,10S,11R)-1,4,4a,4b,5,8,8a,8b-octahydro-1,4-ethano-5,8-methanobiphenylene-10,11-dicarboxylate",
            "smiles": "O=C(OC)[C@@H]1[C@@H](C=C2)[C@@H]3[C@@H]([C@H]4[C@@H]3[C@H]5C=C[C@@H]4C5)[C@@H]2[C@@H]1C(OC)=O"
        },
        'B': {
            "name": "dimethyl (1R,4S,4aS,4bR,5S,8R,8aS,8bR,10S,11R)-1,4,4a,4b,5,8,8a,8b-octahydro-1,4-ethano-5,8-methanobiphenylene-10,11-dicarboxylate",
            "smiles": "O=C(OC)[C@H]1[C@@H](C=C2)[C@@H]3[C@@H]([C@@H]4[C@H]3[C@@H]5C=C[C@H]4C5)[C@@H]2[C@H]1C(OC)=O"
        },
        'C': {
            "name": "dimethyl (1R,4S,4aS,4bR,5R,8S,8aS,8bR,10R,11R)-1,4,4a,4b,5,8,8a,8b-octahydro-1,4-ethano-5,8-methanobiphenylene-10,11-dicarboxylate",
            "smiles": "O=C(OC)[C@@H]1[C@@H](C=C2)[C@@H]3[C@@H]([C@@H]4[C@H]3[C@H]5C=C[C@@H]4C5)[C@@H]2[C@H]1C(OC)=O"
        },
        'D': {
            "name": "dimethyl (1R,4S,4aR,4bR,5R,8S,8aS,8bS,10S,11R)-1,4,4a,4b,5,8,8a,8b-octahydro-1,4-ethano-5,8-methanobiphenylene-10,11-dicarboxylate",
            "smiles": "O=C(OC)[C@@H]1[C@H](C=C2)[C@@H]3[C@@H]([C@H]4[C@@H]3[C@@H]5C=C[C@H]4C5)[C@H]2[C@@H]1C(OC)=O"
        }
    }
    llm_final_answer = 'A'

    error_messages = []

    # --- Verification Step 1: Determine cis/trans configuration from IUPAC names ---
    # The reaction starts with maleic anhydride (a cis-dienophile), so the dicarboxylate groups must remain cis.
    # Cis configuration corresponds to (R,S) or (S,R) at C10 and C11. Trans is (R,R) or (S,S).
    def get_ester_config(name):
        try:
            # Find all R/S descriptors in the name string
            descriptors = re.findall(r'(\d+[a-zA-Z]*[RS])', name)
            # Isolate the descriptors for carbons 10 and 11
            c10 = next(d for d in descriptors if d.startswith('10'))[-1]
            c11 = next(d for d in descriptors if d.startswith('11'))[-1]
            # If descriptors are the same (R,R or S,S), it's trans. Otherwise, it's cis.
            return "trans" if c10 == c11 else "cis"
        except (StopIteration, IndexError):
            return "Unknown"

    ester_configs = {opt: get_ester_config(data['name']) for opt, data in options.items()}
    
    # --- Verification Step 2: Check the LLM's reasoning against the facts ---
    
    # Fact-check claim about option D being trans
    if ester_configs.get('D') == 'cis':
        true_trans_isomer = [opt for opt, config in ester_configs.items() if config == 'trans']
        if true_trans_isomer:
            error_messages.append(
                f"The LLM's reasoning is flawed. It states that option D has 'trans' esters ('In structure D, the two ester groups are trans...'), but the IUPAC name for D (...10S,11R...) indicates a 'cis' configuration. "
                f"The actual 'trans' isomer is option {true_trans_isomer[0]} (...10R,11R...)."
            )

    # Fact-check claim about B and C being identical
    try:
        mol_b = Chem.MolFromSmiles(options['B']['smiles'])
        mol_c = Chem.MolFromSmiles(options['C']['smiles'])
        if mol_b and mol_c:
            canon_b = Chem.MolToSmiles(mol_b, isomericSmiles=True, canonical=True)
            canon_c = Chem.MolToSmiles(mol_c, isomericSmiles=True, canonical=True)
            if canon_b != canon_c:
                error_messages.append(
                    "The LLM's reasoning is flawed. It claims that 'The IUPAC names and SMILES provided for B and C are identical...'. However, their IUPAC names, given SMILES strings, and canonical SMILES are all distinct."
                )
    except Exception:
        error_messages.append("Failed to process SMILES for options B and C with RDKit.")

    # Fact-check the LLM's python code output block
    try:
        canon_b = Chem.MolToSmiles(Chem.MolFromSmiles(options['B']['smiles']), True)
        canon_c = Chem.MolToSmiles(Chem.MolFromSmiles(options['C']['smiles']), True)
        canon_d = Chem.MolToSmiles(Chem.MolFromSmiles(options['D']['smiles']), True)
        # The LLM's output shows B, C, and D having the same canonical SMILES. We check if this is true.
        if not (canon_b == canon_c and canon_c == canon_d):
             error_messages.append("The python code block in the LLM's response produces an incorrect output. The canonical SMILES for options B, C, and D are distinct, not identical as shown in the LLM's output block.")
    except Exception:
        error_messages.append("Failed to generate canonical SMILES to verify the LLM's python output.")


    # --- Final Conclusion ---
    # The chemical principles stated by the LLM (cis-esters required, anti-addition favored) are correct.
    # The final answer 'A' is also correct, as it is the cis-ester, anti-addition product.
    # However, the justification provided by the LLM to reach that answer is factually incorrect.
    
    if error_messages:
        # The reasoning is wrong, even if the final letter is correct.
        return f"Incorrect. While the final answer '{llm_final_answer}' is chemically correct, the provided reasoning is flawed and contains multiple factual errors: " + " ".join(error_messages)
    else:
        # This case is unlikely given the errors found.
        # If no errors were found, we check if the final answer is consistent with the primary constraint.
        if ester_configs.get(llm_final_answer) == 'cis':
            return "Correct"
        else:
            return f"Incorrect. The final answer '{llm_final_answer}' is wrong because it has a '{ester_configs.get(llm_final_answer)}' ester configuration, but a 'cis' configuration is required by the reaction mechanism."

# The function is defined. To get the result, we would call it.
# print(check_correctness_of_llm_answer())