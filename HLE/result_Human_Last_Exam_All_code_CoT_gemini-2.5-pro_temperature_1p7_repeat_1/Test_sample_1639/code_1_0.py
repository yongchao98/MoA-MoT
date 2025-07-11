def calculate_modification_mass():
    """
    Calculates the mass of the variable modification on cysteine.
    The overall process involves:
    1. A complex probe, 2-methylene-4-oxo-4-(prop-2-yn-1-ylamino)butanoic acid, reacts with cysteine.
    2. The protein is enriched and then treated with formic acid.
    3. This treatment cleaves the probe, leaving an itaconic acid adduct (C5H6O4) on the cysteine.
    4. This base adduct has a mass of ~130 Da.
    5. A likely subsequent biological reduction (+4H) results in a final adduct of C5H10O4.
    6. This code calculates the mass of this final adduct C5H10O4.
    """
    # Monoisotopic atomic masses
    atomic_masses = {
        'C': 12.000000,
        'H': 1.007825,
        'O': 15.994915,
        'N': 14.003074,
        'S': 31.972071
    }

    # Formula of the final proposed adduct on cysteine after probe reaction,
    # cleavage, and biological reduction (Itaconic acid + 4H)
    formula = {'C': 5, 'H': 10, 'O': 4}

    total_mass = 0
    equation_parts = []

    for atom, count in formula.items():
        mass = atomic_masses[atom]
        term_mass = count * mass
        total_mass += term_mass
        equation_parts.append(f"{count} * {mass:.6f} ({atom})")

    # Print the explanation and calculation
    print("The final modification is proposed to be a reduced form of an itaconic acid adduct (C5H10O4).")
    print("The mass calculation is as follows:")
    
    # Joining the parts of the equation for printing
    equation_str = " + ".join(equation_parts)
    print(f"Mass = {equation_str} = {total_mass:.6f} Da")

    print(f"\nThe calculated mass is approximately {round(total_mass)} Da.")
    print("This corresponds to answer choice C.")

calculate_modification_mass()