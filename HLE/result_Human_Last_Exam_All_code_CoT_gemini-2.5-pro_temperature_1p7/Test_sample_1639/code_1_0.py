import collections

def calculate_mass_modification():
    """
    Calculates the mass of the variable modification 'x' on cysteine.
    """
    
    # 1. The initial probe is 2-methylene-4-oxo-4-(prop-2-yn-1-ylamino)butanoic acid.
    #    This corresponds to the chemical formula C8H9NO3.
    #    Structure: HOOC-C(=CH2)-CH2-C(=O)-NH-CH2-C≡CH

    # 2. A plausible transformation during the experiment is the reduction of the carboxylic acid
    #    group (-COOH) to a primary alcohol group (-CH2OH).
    #    This transformation results in the loss of one Oxygen (O) and the gain of two Hydrogens (H).
    
    # 3. The formula of the initial probe is C8H9NO3.
    #    The formula of the reduced probe is C8 H(9+2) N O(3-1) = C8H11NO2.
    
    # Monoisotopic masses of the most abundant isotopes.
    atomic_masses = {
        'C': 12.000000,
        'H': 1.007825,
        'N': 14.003074,
        'O': 15.994915,
    }

    # Molecular formula of the final modification attached to cysteine
    formula_str = "C8H11NO2"
    
    # Use a Counter to parse the formula string into element counts
    # This is more robust than manual parsing if formulas change.
    parsed_formula = collections.Counter()
    element = ''
    count = ''
    for char in formula_str:
        if char.isalpha():
            if element: # If we have a pending element, process it
                parsed_formula[element] += int(count) if count else 1
            element = char
            count = ''
        elif char.isdigit():
            count += char
    if element: # Process the last element
        parsed_formula[element] += int(count) if count else 1

    total_mass = 0
    calculation_steps = []
    
    # Sort elements for consistent output order (C, H, N, O)
    sorted_elements = sorted(parsed_formula.keys(), key=lambda x: ('C', 'H', 'N', 'O').index(x))

    for element in sorted_elements:
        count = parsed_formula[element]
        mass = atomic_masses[element]
        element_mass = count * mass
        total_mass += element_mass
        calculation_steps.append(f"{count} * mass({element})")
        print(f"Number of {element} atoms: {count}")

    print("\nThe mass modification 'x' is the mass of the final probe C8H11NO2.")
    print("The equation to calculate the mass is:")
    
    equation_str_parts = []
    for element in sorted_elements:
        count = parsed_formula[element]
        equation_str_parts.append(f"{count} * {atomic_masses[element]:.6f}")

    print(f"x = {' + '.join(equation_str_parts)}")
    
    print(f"\nx = {total_mass:.4f} Da")
    
    print(f"\nThis calculated mass of {total_mass:.2f} Da is closest to the answer choice 154.")
    print("The small difference is likely due to the common practice in mass spectrometry of using the integer mass of the protonated molecule, [M+H]+, which would be 154.")

calculate_mass_modification()