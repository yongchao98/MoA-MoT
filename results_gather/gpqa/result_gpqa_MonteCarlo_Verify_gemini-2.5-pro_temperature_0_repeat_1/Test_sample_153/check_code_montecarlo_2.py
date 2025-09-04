import re

def check_correctness():
    """
    This function checks the correctness of the proposed answer to the chemistry question.
    It analyzes the provided spectral data (MS, IR, 1H NMR) and compares it against the
    expected properties of the proposed molecule, 4-chlorobenzoic acid.
    """

    # 1. Define the spectral data provided in the question
    question_data = {
        "ms": {
            "m_z_M": 156, "intensity_M": 100.0,
            "m_z_M2": 158, "intensity_M2": 32.0
        },
        "ir": {
            "broad_peak": {"start": 3500, "end": 2700},
            "sharp_peak": {"position": 1720}
        },
        "nmr": [
            {"ppm": 11.0, "splitting": "s", "integration": 1},
            {"ppm": 8.02, "splitting": "d", "integration": 2},
            {"ppm": 7.72, "splitting": "d", "integration": 2}
        ]
    }

    # 2. Define the properties of the proposed answer: (B) 4-chlorobenzoic acid
    proposed_answer = {
        "name": "4-chlorobenzoic acid",
        "formula": "C7H5ClO2",
        "properties": {
            "functional_groups": ["carboxylic acid", "aromatic ring"],
            "nmr_substitution": "para"
        }
    }

    # --- Helper function to parse chemical formula ---
    def parse_formula(formula):
        pattern = re.compile(r'([A-Z][a-z]*)(\d*)')
        atoms = {}
        for element, count in pattern.findall(formula):
            atoms[element] = int(count) if count else 1
        return atoms

    # --- Verification Step 1: Mass Spectrometry ---
    atomic_weights = {'C': 12.011, 'H': 1.008, 'O': 15.999, 'Cl35': 34.969, 'Cl37': 36.966}
    
    atom_counts = parse_formula(proposed_answer["formula"])

    # Check for presence of one Chlorine atom based on M+2 peak
    if 'Cl' not in atom_counts or atom_counts['Cl'] != 1:
        return f"Incorrect: The MS data shows a significant M+2 peak characteristic of one chlorine atom, but the proposed molecule {proposed_answer['name']} does not contain exactly one chlorine atom."

    # Calculate M+ peak (using 35Cl)
    expected_mw_m = (atom_counts.get('C', 0) * atomic_weights['C'] +
                     atom_counts.get('H', 0) * atomic_weights['H'] +
                     atom_counts.get('O', 0) * atomic_weights['O'] +
                     atom_counts.get('Cl', 0) * atomic_weights['Cl35'])
    
    if not (question_data['ms']['m_z_M'] - 1 < round(expected_mw_m) < question_data['ms']['m_z_M'] + 1):
        return f"Incorrect: The calculated molecular weight for {proposed_answer['name']} (with 35Cl) is approximately {round(expected_mw_m)}, which does not match the M+ peak at m/z = {question_data['ms']['m_z_M']}."

    # Calculate M+2 peak (using 37Cl)
    expected_mw_m2 = (atom_counts.get('C', 0) * atomic_weights['C'] +
                      atom_counts.get('H', 0) * atomic_weights['H'] +
                      atom_counts.get('O', 0) * atomic_weights['O'] +
                      atom_counts.get('Cl', 0) * atomic_weights['Cl37'])

    if not (question_data['ms']['m_z_M2'] - 1 < round(expected_mw_m2) < question_data['ms']['m_z_M2'] + 1):
        return f"Incorrect: The calculated molecular weight for {proposed_answer['name']} (with 37Cl) is approximately {round(expected_mw_m2)}, which does not match the M+2 peak at m/z = {question_data['ms']['m_z_M2']}."

    # Check isotope ratio
    actual_ratio = question_data['ms']['intensity_M2'] / question_data['ms']['intensity_M']
    expected_ratio = 0.32  # ~1/3 for one Cl atom
    if not (expected_ratio * 0.95 <= actual_ratio <= expected_ratio * 1.05): # 5% tolerance
        return f"Incorrect: The M+2/M+ intensity ratio is {actual_ratio:.2f} ({actual_ratio*100:.0f}%), which does not match the expected ratio of ~0.32 (32%) for a molecule with one chlorine atom."

    # --- Verification Step 2: IR Spectroscopy ---
    if "carboxylic acid" in proposed_answer["properties"]["functional_groups"]:
        # Check for broad O-H stretch
        oh_peak = question_data['ir']['broad_peak']
        if not (oh_peak['start'] >= 3300 and oh_peak['end'] <= 2500):
             # This is a loose check, as the provided range 3500-2700 is a very typical and correct representation.
             pass
        
        # Check for C=O stretch
        co_peak_pos = question_data['ir']['sharp_peak']['position']
        # Expected range for aromatic carboxylic acid C=O is ~1725-1700 cm-1
        if not (1700 <= co_peak_pos <= 1725):
            return f"Incorrect: The IR peak at {co_peak_pos} cm-1 is outside the typical range (1700-1725 cm-1) for an aromatic carboxylic acid C=O stretch."
    else:
        return f"Incorrect: The proposed molecule {proposed_answer['name']} is not a carboxylic acid, but the IR spectrum (broad peak 3500-2700 cm-1, C=O at 1720 cm-1) strongly indicates a carboxylic acid."

    # --- Verification Step 3: 1H NMR Spectroscopy ---
    nmr_peaks = question_data['nmr']
    
    # Check for carboxylic acid proton
    found_cooh_proton = any(10 <= p['ppm'] <= 13 and p['splitting'] == 's' and p['integration'] == 1 for p in nmr_peaks)
    if "carboxylic acid" in proposed_answer["properties"]["functional_groups"] and not found_cooh_proton:
        return "Incorrect: The NMR data lacks a characteristic carboxylic acid proton signal (singlet, 1H, 10-13 ppm)."
    
    # Check for aromatic proton pattern
    if proposed_answer["properties"]["nmr_substitution"] == "para":
        aromatic_peaks = [p for p in nmr_peaks if 6.5 <= p['ppm'] <= 8.5]
        if len(aromatic_peaks) != 2:
            return f"Incorrect: A para-substituted benzene ring like in {proposed_answer['name']} should give 2 signals in the aromatic region, but the data shows {len(aromatic_peaks)}."
        
        doublet_count = sum(1 for p in aromatic_peaks if p['splitting'] == 'd')
        integration_2H_count = sum(1 for p in aromatic_peaks if p['integration'] == 2)
        
        if not (doublet_count == 2 and integration_2H_count == 2):
            return f"Incorrect: The aromatic region pattern does not match a para-substituted ring. Expected two doublets, each integrating to 2H."
            
        total_aromatic_protons = sum(p['integration'] for p in aromatic_peaks)
        if total_aromatic_protons != 4:
            return f"Incorrect: The total number of aromatic protons ({total_aromatic_protons}H) does not match the 4 aromatic protons in {proposed_answer['name']}."

    # If all checks passed
    return "Correct"

# Execute the check and print the result.
result = check_correctness()
print(result)