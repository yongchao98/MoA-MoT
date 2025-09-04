import sys
from io import StringIO

def check_chemistry_answer():
    """
    This function checks the correctness of the final answer for the given chemistry problem.
    It analyzes the spectral data step-by-step to determine the correct structure and then
    compares it with the provided answer.
    """
    
    # --- 1. Define Candidate Structures based on the Question's Options ---
    # This mapping is crucial and must match the question exactly.
    candidates = {
        'A': {
            'name': '4-chlorobenzoic acid',
            'mw_35cl': 156,
            'has_chlorine': True,
            'functional_group': 'carboxylic_acid',
            'substitution': 'para',
            'nmr_aromatic_protons': 4,
            'nmr_aromatic_pattern_desc': "Symmetrical 1,4-(para) substitution gives two doublets, each for 2H."
        },
        'B': {
            'name': '3-Chloro-2-hydroxybenzaldehyde',
            'mw_35cl': 156,
            'has_chlorine': True,
            'functional_group': 'aldehyde_phenol',
            'substitution': '1,2,3-trisubstituted',
            'nmr_aromatic_protons': 3,
            'nmr_aromatic_pattern_desc': "Has only 3 aromatic protons, not 4."
        },
        'C': {
            'name': '2-chlorobenzoic acid',
            'mw_35cl': 156,
            'has_chlorine': True,
            'functional_group': 'carboxylic_acid',
            'substitution': 'ortho',
            'nmr_aromatic_protons': 4,
            'nmr_aromatic_pattern_desc': "Unsymmetrical 1,2-(ortho) substitution would give a complex pattern for 4 unique protons, not two simple doublets."
        },
        'D': {
            'name': 'Phenyl chloroformate',
            'mw_35cl': 156,
            'has_chlorine': True,
            'functional_group': 'chloroformate',
            'substitution': 'monosubstituted',
            'nmr_aromatic_protons': 5,
            'nmr_aromatic_pattern_desc': "Has 5 aromatic protons, not 4, and no carboxylic acid proton."
        }
    }

    # --- 2. Define Spectral Data Constraints from the Question ---
    constraints = {
        'ms_m_plus_2_ratio': 0.32,  # ~3:1 ratio for one Cl atom
        'ir_carboxylic_acid': True, # Broad 3500-2700 cm-1 and sharp 1720 cm-1
        'nmr_acid_proton': True,    # 11.0 ppm (s, 1H)
        'nmr_aromatic_pattern': '2 doublets, 2H each' # 8.02 ppm (d, 2H), 7.72 (d, 2H)
    }

    # --- 3. Systematically Determine the Correct Answer ---
    survivors = list(candidates.keys())
    reasons_for_elimination = {}

    # Constraint: Mass Spectrometry (All candidates pass this)
    # The M+2 peak at 32% confirms one chlorine atom, which all candidates have.
    # The MW of 156 is also consistent with all candidates.
    
    # Constraint: IR Spectrum (identifies carboxylic acid)
    survivors_after_ir = []
    for option in survivors:
        if candidates[option]['functional_group'] == 'carboxylic_acid':
            survivors_after_ir.append(option)
        else:
            reasons_for_elimination[option] = f"Fails IR constraint. It is a {candidates[option]['functional_group'].replace('_', ' ')}, not a carboxylic acid as indicated by the broad 3500-2700 cm-1 and sharp 1720 cm-1 peaks."
    survivors = survivors_after_ir

    # Constraint: NMR Spectrum (identifies substitution pattern)
    survivors_after_nmr = []
    for option in survivors:
        # Check for carboxylic acid proton (all remaining candidates have this)
        # Check aromatic pattern
        if candidates[option]['substitution'] == 'para':
             survivors_after_nmr.append(option)
        else:
            reasons_for_elimination[option] = f"Fails NMR constraint. {candidates[option]['nmr_aromatic_pattern_desc']}"
    survivors = survivors_after_nmr
    
    # --- 4. Final Verification ---
    # The provided answer to check
    llm_answer = 'A'

    if not survivors:
        return "Error in analysis: No candidate satisfies all constraints."
    
    ground_truth_answer = survivors[0]

    if llm_answer == ground_truth_answer:
        return "Correct"
    else:
        # Build a detailed error message
        output = StringIO()
        print(f"Incorrect. The provided answer was <<<{llm_answer}>>>, but the correct answer is <<<{ground_truth_answer}>>>.", file=output)
        
        llm_candidate_info = candidates[llm_answer]
        correct_candidate_info = candidates[ground_truth_answer]

        print("\nReasoning:", file=output)
        print(f"1. The correct structure is {correct_candidate_info['name']} ({ground_truth_answer}) because it is the only one that satisfies all spectral data.", file=output)
        print(f"   - MS: Has MW of 156 and one Cl atom.", file=output)
        print(f"   - IR: Is a carboxylic acid, matching the broad O-H and C=O peaks.", file=output)
        print(f"   - NMR: Is para-substituted, which perfectly explains the symmetrical pattern of two doublets for 2H each in the aromatic region.", file=output)
        
        print(f"\n2. The proposed answer {llm_candidate_info['name']} ({llm_answer}) is incorrect.", file=output)
        if llm_answer in reasons_for_elimination:
            print(f"   - It fails the following constraint: {reasons_for_elimination[llm_answer]}", file=output)
        else:
            # This case shouldn't happen with this problem, but is good practice
            print("   - It was eliminated during the analysis.", file=output)

        return output.getvalue()

# Run the check
result = check_chemistry_answer()
print(result)