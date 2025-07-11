import sys

def solve_helix_type():
    """
    This function determines the likely helix type for an alternating alpha/epsilon-peptide
    by extrapolating from known data for alpha/beta and alpha/gamma peptides.
    """

    # Step 1 & 2: Define backbone sizes for the amino acid types.
    # The backbone is defined as the chain of atoms from the amide N to the carbonyl C'.
    # Alpha-amino acid (e.g., Alanine): -NH-C(alpha)-CO-
    s_alpha = 3
    # Beta-amino acid: -NH-C(beta)-C(alpha)-CO-
    s_beta = 4
    # Gamma-amino acid: -NH-C(gamma)-C(beta)-C(alpha)-CO-
    s_gamma = 5
    # Epsilon-amino acid: -NH-(CH2)5-CO-
    s_epsilon = 7

    print("Step 1: Establishing the properties of known foldamer helices.")
    # Step 3: Define known helix types for alternating copolymers.
    # For an alpha/X peptide with i->i+2 H-bonds, we have two ring types:
    # m_X: The ring formed by an X->X H-bond, which spans an alpha residue.
    # m_alpha: The ring formed by an alpha->alpha H-bond, which spans an X residue.
    
    # Data point 1: alpha/beta-peptide
    ab_peptide_data = {'s_X': s_beta, 'm_X': 11, 'm_alpha': 9}
    print(f"  - For an alpha/beta-peptide (backbone sizes {s_alpha}, {s_beta}), the helix is an {ab_peptide_data['m_X']}/{ab_peptide_data['m_alpha']}-helix.")

    # Data point 2: alpha/gamma-peptide
    ag_peptide_data = {'s_X': s_gamma, 'm_X': 14, 'm_alpha': 12}
    print(f"  - For an alpha/gamma-peptide (backbone sizes {s_alpha}, {s_gamma}), the helix is a {ag_peptide_data['m_X']}/{ag_peptide_data['m_alpha']}-helix.")
    
    # Python 2 compatibility for printing
    if sys.version_info[0] < 3:
      print ''

    print("\nStep 2: Finding the linear relationship.")
    # Step 4: Calculate the rate of change (slope) for the ring sizes
    # as the backbone size of the 'X' residue increases.
    delta_s_X = ag_peptide_data['s_X'] - ab_peptide_data['s_X']
    delta_m_X = ag_peptide_data['m_X'] - ab_peptide_data['m_X']
    delta_m_alpha = ag_peptide_data['m_alpha'] - ab_peptide_data['m_alpha']

    slope_m_X = delta_m_X / float(delta_s_X)
    slope_m_alpha = delta_m_alpha / float(delta_s_X)
    
    print(f"When the backbone size of the non-alpha residue increases by {delta_s_X} (from beta to gamma),")
    print(f"  - The ring spanning the alpha-residue (m_X) increases by {delta_m_X}.")
    print(f"  - The ring spanning the non-alpha residue (m_alpha) increases by {delta_m_alpha}.")
    print(f"This indicates a consistent slope of {int(slope_m_alpha)} atoms per unit increase in backbone size.")
    if sys.version_info[0] < 3:
      print ''


    print("\nStep 3: Predicting the helix type for the alpha/epsilon-peptide.")
    # Step 5: Extrapolate to the alpha/epsilon-peptide case.
    target_s_X = s_epsilon
    print(f"The backbone size of the epsilon-amino acid is {target_s_X}.")

    # Use the point-slope form: y - y1 = m * (x - x1) => y = m * (x - x1) + y1
    # We use the alpha/beta peptide as our base point (x1, y1).
    
    # Predict the size of the ring spanning the alpha-residue (m_epsilon)
    predicted_m_X = int(slope_m_X * (target_s_X - ab_peptide_data['s_X']) + ab_peptide_data['m_X'])
    print("Predicted ring size (m_epsilon, spanning an alpha-residue):")
    print(f"  m_epsilon = {int(slope_m_X)} * ({target_s_X} - {ab_peptide_data['s_X']}) + {ab_peptide_data['m_X']} = {predicted_m_X}")
    
    # Predict the size of the ring spanning the epsilon-residue (m_alpha)
    predicted_m_alpha = int(slope_m_alpha * (target_s_X - ab_peptide_data['s_X']) + ab_peptide_data['m_alpha'])
    print("Predicted ring size (m_alpha, spanning an epsilon-residue):")
    print(f"  m_alpha = {int(slope_m_alpha)} * ({target_s_X} - {ab_peptide_data['s_X']}) + {ab_peptide_data['m_alpha']} = {predicted_m_alpha}")
    
    if sys.version_info[0] < 3:
      print ''
    print(f"\nThe predicted helix type is therefore a {predicted_m_alpha}/{predicted_m_X}-helix.")


solve_helix_type()