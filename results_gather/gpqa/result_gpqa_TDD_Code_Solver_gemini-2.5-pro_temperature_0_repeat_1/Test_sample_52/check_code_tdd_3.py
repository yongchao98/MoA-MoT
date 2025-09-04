import re

def check_correctness():
    """
    This function checks the correctness of the provided answer by systematically verifying
    the constraints given in the chemistry question.

    The question describes a di-substituted 6-membered aromatic ring compound with the following features:
    1.  FTIR: Ester group present.
    2.  1H NMR:
        - 2 aromatic-H signals
        - 2 vinyl-H signals (one doublet, one doublet of quartets)
        - 2 -CH3 group signals
        - 0 -CH2 group signals
    3.  The proposed answer is option A, which corresponds to the formula C11H12O2.
    """
    
    # The molecular formula from the selected answer 'A'
    formula_to_check = "C11H12O2"

    # --- Step 1: Parse the molecular formula into atom counts ---
    try:
        c_match = re.search(r'C(\d+)', formula_to_check)
        h_match = re.search(r'H(\d+)', formula_to_check)
        o_match = re.search(r'O(\d+)', formula_to_check)
        
        C = int(c_match.group(1))
        H = int(h_match.group(1))
        O = int(o_match.group(1))
    except (AttributeError, ValueError):
        return f"Error: Could not parse the chemical formula '{formula_to_check}'."

    # --- Step 2: Deduce structural fragments from the spectral data ---
    
    # Fragment 1: A di-substituted 6-membered aromatic ring has a core of C6H4.
    # This accounts for the aromatic-H signals.
    fragments = {'ring': {'C': 6, 'H': 4, 'O': 0}}
    
    # Fragment 2: An ester group (-COO-) is present (from FTIR).
    # This requires at least 2 oxygen atoms.
    if O < 2:
        return f"Incorrect. The formula {formula_to_check} has {O} oxygen(s), but an ester group requires at least 2."

    # Fragment 3: NMR data analysis
    # - The two vinyl-H signals (one doublet, one doublet of quartets) combined with one of the -CH3 signals
    #   is a classic pattern for a propenyl group (-CH=CH-CH3).
    # - The absence of -CH2 signals supports this and rules out groups like ethyl or propyl.
    # - The second -CH3 signal must belong to the other part of the ester. The simplest arrangement is a methyl ester,
    #   where the ester group is a methyl carboxylate (-COOCH3).
    
    # Let's define the two substituents based on this deduction:
    # Substituent 1: Propenyl group -> -C3H5
    fragments['substituent_1'] = {'C': 3, 'H': 5, 'O': 0}
    
    # Substituent 2: Methyl carboxylate group -> -COOCH3
    fragments['substituent_2'] = {'C': 2, 'H': 3, 'O': 2}

    # --- Step 3: Assemble the fragments and verify the total formula ---
    
    c_calculated = sum(f['C'] for f in fragments.values())
    h_calculated = sum(f['H'] for f in fragments.values())
    o_calculated = sum(f['O'] for f in fragments.values())
    
    # Check if the calculated formula matches the formula from the chosen answer
    if C == c_calculated and H == h_calculated and O == o_calculated:
        # The formula is consistent. As a final check, we can verify the Degrees of Unsaturation (DBU).
        # DBU from formula = C - H/2 + 1
        dbu_from_formula = C - (H / 2) + 1
        
        # DBU from deduced structure: 4 (aromatic ring) + 1 (propenyl C=C) + 1 (ester C=O) = 6
        expected_dbu = 6
        
        if dbu_from_formula == expected_dbu:
            return "Correct"
        else:
            # This case is unlikely if the atom counts match, but serves as a robust cross-check.
            return f"Incorrect. The formula {formula_to_check} has a DBU of {dbu_from_formula}, but the deduced structure requires a DBU of {expected_dbu}."
    else:
        return f"Incorrect. Assembling the fragments deduced from the spectral data (C6H4 ring, C3H5 propenyl, C2H3O2 methyl ester) results in a formula of C{c_calculated}H{h_calculated}O{o_calculated}, which does not match the proposed formula {formula_to_check}."

# Execute the check
result = check_correctness()
print(result)