def calculate_modification_mass():
    """
    Calculates the mass of the variable modification on cysteine.
    
    1. The probe is an itaconic acid derivative. Itaconic acid is C5H6O4.
    2. It reacts with cysteine via Michael addition, which saturates the double bond.
       This adds two hydrogen atoms to the itaconic acid core.
    3. Formic acid cleaves the probe, leaving a methylsuccinic acid group attached.
    4. The formula of this final modification is C5H8O4.
    5. We calculate the monoisotopic mass of this C5H8O4 group.
    """

    # Monoisotopic masses of the most common isotopes
    atomic_mass = {
        'C': 12.000000,
        'H': 1.007825,
        'O': 15.994915,
    }

    # Formula of the final modification (methylsuccinic acid)
    formula = {'C': 5, 'H': 8, 'O': 4}

    # Calculate the mass
    total_mass = 0
    
    print("Calculating the mass of the final adducted molecule (C5H8O4):")
    
    # C
    c_mass = formula['C'] * atomic_mass['C']
    total_mass += c_mass
    print(f"Mass from {formula['C']} Carbon atoms: {formula['C']} * {atomic_mass['C']} = {c_mass:.4f} Da")

    # H
    h_mass = formula['H'] * atomic_mass['H']
    total_mass += h_mass
    print(f"Mass from {formula['H']} Hydrogen atoms: {formula['H']} * {atomic_mass['H']} = {h_mass:.4f} Da")
    
    # O
    o_mass = formula['O'] * atomic_mass['O']
    total_mass += o_mass
    print(f"Mass from {formula['O']} Oxygen atoms: {formula['O']} * {atomic_mass['O']} = {o_mass:.4f} Da")

    print("\nFinal Equation:")
    print(f"({formula['C']} * {atomic_mass['C']}) + ({formula['H']} * {atomic_mass['H']}) + ({formula['O']} * {atomic_mass['O']}) = {total_mass:.4f}")

    print(f"\nThe calculated mass of the modification is approximately {round(total_mass)} Da.")
    print("Comparing this to the answer choices, the closest value is 134.")

calculate_modification_mass()