import pandas as pd

def propose_methanol_charges():
    """
    This function proposes a reasonable set of partial charges for a new
    all-atom model of methanol (CH3OH) and verifies its neutrality.

    The selection is based on established principles of molecular modeling:
    1. The net charge of the molecule must be zero.
    2. Charges should reflect atomic electronegativity (O < C < H).
    3. Chemically equivalent atoms (the 3 methyl hydrogens) should have identical charges.

    Based on these principles, Choice D is the most reasonable model.
    """

    charges = {
        'Atom': ['Carbon', 'Methyl H1', 'Methyl H2', 'Methyl H3', 'Oxygen', 'Hydroxyl hydrogen'],
        'Symbol': ['C', 'H', 'H', 'H', 'O', 'H'],
        'Partial Charge (e)': [0.1450, 0.0400, 0.0400, 0.0400, -0.6830, 0.4180]
    }

    df = pd.DataFrame(charges)
    total_charge = df['Partial Charge (e)'].sum()

    print("Proposed Partial Charges for Methanol (CH3OH):")
    # Using pandas to_string() to format output without the index
    print(df.to_string(index=False))

    print("\nVerifying charge neutrality...")
    equation_parts = [f"{charge:+.4f}" for charge in df['Partial Charge (e)']]
    equation = " + ".join(equation_parts).replace('+ -', '- ')
    print(f"Sum of charges: {equation} = {total_charge:.4f}")

    if abs(total_charge) < 1e-9:
        print("\nThe molecule is charge-neutral. This is a valid set of charges.")
    else:
        print(f"\nWarning: The molecule is not charge-neutral (Total = {total_charge:.4f}).")

propose_methanol_charges()