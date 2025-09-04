import math

def check_spectroscopy_answer():
    """
    This function programmatically checks if the provided spectral data uniquely identifies
    4-chlorobenzoic acid among the given options. It evaluates the Mass Spec, IR, and NMR data
    against the chemical properties of each candidate molecule.
    """

    # --- Define Spectral Data from the Question ---
    # Mass Spec: m/z = 156 (M+), m/z = 158 (M+2, 32% intensity of M+)
    # IR: Broad peak 3500-2700 cm-1, strong sharp peak 1720 cm-1
    # 1H NMR: 11.0 ppm (s, 1H), 8.02 ppm (d, 2H), 7.72 ppm (d, 2H)

    # --- Define Properties of Candidate Molecules ---
    candidates = {
        'A': {'name': '4-chlorobenzoic acid', 'formula': 'C7H5ClO2', 'functional_group': 'carboxylic_acid', 'substitution': 'para', 'aromatic_protons': 4},
        'B': {'name': '2-chlorobenzoic acid', 'formula': 'C7H5ClO2', 'functional_group': 'carboxylic_acid', 'substitution': 'ortho', 'aromatic_protons': 4},
        'C': {'name': 'Phenyl chloroformate', 'formula': 'C7H5ClO2', 'functional_group': 'ester/acid_chloride', 'substitution': 'mono', 'aromatic_protons': 5},
        'D': {'name': '3-Chloro-2-hydroxybenzaldehyde', 'formula': 'C7H5ClO2', 'functional_group': 'phenol/aldehyde', 'substitution': 'tri-substituted', 'aromatic_protons': 3}
    }
    
    proposed_answer_key = 'A'
    proposed_answer_data = candidates[proposed_answer_key]

    # --- Verification Step 1: Mass Spectrometry ---
    # All candidates have the formula C7H5ClO2. Let's verify the mass.
    # Using most abundant isotopes: C=12, H=1, O=16, Cl=35
    m_plus_calc = 7 * 12 + 5 * 1 + 1 * 35 + 2 * 16
    if m_plus_calc != 156:
        return f"Mass Spec check failed: Calculated M+ peak for C7H5(35)ClO2 is {m_plus_calc}, but the data shows 156."
    
    # Using the heavier chlorine isotope: Cl=37
    m_plus_2_calc = 7 * 12 + 5 * 1 + 1 * 37 + 2 * 16
    if m_plus_2_calc != 158:
        return f"Mass Spec check failed: Calculated M+2 peak for C7H5(37)ClO2 is {m_plus_2_calc}, but the data shows 158."

    # Check the M+2 isotope ratio. Natural abundance of 35Cl is ~75.8% and 37Cl is ~24.2%.
    # The expected ratio of M+2 to M+ is ~24.2 / 75.8 ≈ 0.319.
    # The observed ratio is 32% / 100% = 0.32. This is an excellent match for one chlorine atom.
    if not math.isclose(0.32, 24.2/75.8, abs_tol=0.02):
         return "Mass Spec check failed: The observed M+2 ratio of 32% is not consistent with the presence of one chlorine atom."
    # Conclusion: MS confirms the molecular formula C7H5ClO2, which is consistent with all four options.

    # --- Verification Step 2: IR Spectroscopy ---
    # A very broad peak from 3500-2700 cm-1 combined with a strong C=O peak (~1720 cm-1) is a definitive sign of a carboxylic acid.
    if proposed_answer_data['functional_group'] != 'carboxylic_acid':
        return f"IR check failed: The IR data strongly indicates a carboxylic acid, but the proposed answer, {proposed_answer_data['name']}, is not one."
    
    # The IR data correctly rules out options C (ester/acid chloride) and D (phenol/aldehyde), which do not have this specific combination of peaks.
    if candidates['C']['functional_group'] == 'carboxylic_acid' or candidates['D']['functional_group'] == 'carboxylic_acid':
        return "Constraint check failed: The logic that IR rules out options C and D is flawed because the checker's internal data is incorrect."
    # Conclusion: IR narrows the possibilities to A and B.

    # --- Verification Step 3: 1H NMR Spectroscopy ---
    # Signal 1: 11.0 ppm (s, 1H). This is the acidic proton of a carboxylic acid. Consistent with A and B.
    # Aromatic signals: 8.02 ppm (d, 2H) and 7.72 ppm (d, 2H).
    # This pattern of two doublets, each integrating to two protons, is characteristic of a symmetrical, 1,4-disubstituted (para) benzene ring.
    
    # Check if the proposed answer's structure matches this pattern.
    if proposed_answer_data['substitution'] != 'para':
        return f"NMR check failed: The NMR aromatic pattern indicates para-substitution, but the proposed answer, {proposed_answer_data['name']}, has a '{proposed_answer_data['substitution']}' substitution."

    # Check if the other remaining option (B) is correctly ruled out.
    # 2-chlorobenzoic acid is ortho-substituted and unsymmetrical. It would show four distinct 1H signals in the aromatic region, not two 2H signals.
    if candidates['B']['substitution'] == 'para':
        return "Constraint check failed: The logic that NMR rules out option B is flawed. 2-chlorobenzoic acid is ortho-substituted, not para."
        
    # Final check: All pieces of evidence point to 4-chlorobenzoic acid.
    # MS confirms the formula. IR confirms the carboxylic acid functional group. NMR confirms the para-substitution pattern.
    # The provided explanation correctly uses this logic to arrive at the answer.
    
    return "Correct"

# Run the check
result = check_spectroscopy_answer()
print(result)