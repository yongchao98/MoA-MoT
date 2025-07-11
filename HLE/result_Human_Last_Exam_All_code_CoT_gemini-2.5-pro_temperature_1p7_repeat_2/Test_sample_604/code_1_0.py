import pandas as pd

def analyze_hyperfine_field():
    """
    Analyzes different Fe species to find which one is expected to have the largest hyperfine field.
    The hyperfine field B_hf is modeled by the equation:
    B_hf = B_contact + B_orbital
    where B_contact is proportional to the spin S and B_orbital reduces the field if L is not quenched.
    A positive value here denotes the magnitude of the field.
    """

    # Data for the given choices.
    # L_quenched is True for orbitally non-degenerate ground states (like high-spin d5).
    choices = [
        {'option': 'A', 'species': 'sq. pyramidal Fe(II)', 'spin': 0, 'L_quenched': True},
        {'option': 'B', 'species': 'planar Fe(III)', 'spin': 2.5, 'L_quenched': True}, # High-spin d5 has a non-degenerate 6A1 ground state, so L is quenched.
        {'option': 'C', 'species': 'linear Fe(II)', 'spin': 2, 'L_quenched': False}, # High-spin d6 in linear geometry is orbitally degenerate (e.g., 5-Delta_g)
        {'option': 'D', 'species': 'tetrahedral Fe(II)', 'spin': 2, 'L_quenched': False},# High-spin d6 in Td has an orbitally degenerate 5E ground state.
        {'option': 'E', 'species': 'trig. bipy. Fe(IV)', 'spin': 2, 'L_quenched': False} # High-spin d4 in D3h has an orbitally degenerate 5E' ground state.
    ]

    # Approximate physical constants for an illustrative calculation (in Tesla).
    # The Fermi contact field per unit of spin for Fe is roughly 22 T.
    # The orbital field contribution is complex, but when unquenched it opposes the contact term.
    # We use an illustrative reduction value.
    CONTACT_TERM_PER_SPIN = 22.0  # T per spin unit
    ORBITAL_REDUCTION = 15.0       # T (illustrative value when L is unquenched)

    print("Analysis of Hyperfine Field (B_hf) for each option:")
    print("-----------------------------------------------------")
    print("Formula: |B_hf| ~ |(S * C_contact) - B_orbital|")
    print(f"Constants used: C_contact = {CONTACT_TERM_PER_SPIN} T/spin, B_orbital = {ORBITAL_REDUCTION} T if L is unquenched\n")

    results = []
    for choice in choices:
        s_val = choice['spin']
        is_l_quenched = choice['L_quenched']

        # Calculate components of the hyperfine field
        b_contact = s_val * CONTACT_TERM_PER_SPIN
        b_orbital = 0.0 if is_l_quenched else ORBITAL_REDUCTION

        # The orbital field opposes the contact field
        b_hf_magnitude = abs(b_contact - b_orbital)

        result = {
            "Option": choice['option'],
            "Species": choice['species'],
            "S": s_val,
            "L Quenched?": "Yes" if is_l_quenched else "No",
            "B_contact (T)": b_contact,
            "B_orbital (T)": b_orbital,
            "|B_hf| (T)": b_hf_magnitude
        }
        results.append(result)

    # Use pandas for pretty printing the table
    df = pd.DataFrame(results)
    print(df.to_string(index=False))

    # Determine the best option
    best_option = max(results, key=lambda x: x['|B_hf| (T)'])

    print("\nConclusion:")
    print(f"Option {best_option['Option']} is expected to have the largest hyperfine field.")
    print("This is because it combines the maximum possible spin state (S = 5/2) with a quenched orbital angular momentum (L=0),")
    print("maximizing the dominant Fermi contact contribution without reduction from an opposing orbital field.")
    print("\nIllustrative final equation for the winning option:")
    winning_b_contact = best_option['S'] * CONTACT_TERM_PER_SPIN
    winning_b_orbital = 0
    print(f"|B_hf({best_option['Option']})| ~ |({best_option['S']} * {CONTACT_TERM_PER_SPIN}) - {winning_b_orbital}| = {best_option['|B_hf| (T)|']} T")


if __name__ == "__main__":
    analyze_hyperfine_field()
<<<B>>>