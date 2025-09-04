import collections

def check_answer():
    """
    Checks the correctness of the answer to the Diels-Alder NOESY problem.
    """
    # 1. Define the proton signals from the options
    signals = {
        "bridgehead_Me": {"integration": "6H", "multiplicity": "singlet", "shift_ppm": "~1"},
        "vinylic_Me": {"integration": "6H", "multiplicity": "singlet", "shift_ppm": "~1.7"},
        "bridge_H": {"integration": "1H", "multiplicity": "doublet", "shift_ppm": "~1.5"},
        "anhydride_H": {"integration": "2H", "multiplicity": "singlet", "shift_ppm": "~3.5"},
    }

    # 2. Define the options based on the question text
    options = {
        "A": {"bridgehead_Me", "vinylic_Me"},
        "B": {"bridge_H", "anhydride_H"},
        "C": {"vinylic_Me", "anhydride_H"},
        "D": {"bridgehead_Me", "bridge_H"},
    }
    
    final_answer_from_llm = "C"

    # 3. Model the chemical principles
    # Principle 1: The diene (1,2,3,4-tetramethyl-1,3-cyclopentadiene) is highly sterically hindered.
    is_sterically_hindered = True

    # Principle 2: The Alder-endo rule is overridden by severe steric hindrance.
    # The exo product becomes the major product.
    if is_sterically_hindered:
        major_product_type = "exo"
        minor_product_type = "endo"
    else:
        # Default endo rule for non-hindered systems
        major_product_type = "endo"
        minor_product_type = "exo"

    # Principle 3: Define spatial proximities (<5 Angstroms) for NOESY cross-peaks in each isomer.
    # This is a simplified model of the 3D structures.
    # In the exo isomer, anhydride protons are close to vinylic methyls.
    exo_proximities = {("anhydride_H", "vinylic_Me"), ("bridgehead_Me", "bridge_H")}
    # In the endo isomer, anhydride protons are close to the C7 bridge protons.
    endo_proximities = {("anhydride_H", "bridge_H"), ("bridgehead_Me", "bridge_H")}

    # 4. Determine the unique cross-peak
    # The question asks for a cross-peak present in the major product but absent/weak in the minor.
    if major_product_type == "exo":
        major_peaks = exo_proximities
        minor_peaks = endo_proximities
    else: # major_product_type == "endo"
        major_peaks = endo_proximities
        minor_peaks = exo_proximities

    # Find the set of peaks that are in the major product but not the minor
    unique_major_peaks = major_peaks - minor_peaks

    if not unique_major_peaks:
        return "Logic Error: Could not find a unique cross-peak distinguishing the major and minor products based on the model."

    # The problem implies a single distinguishing cross-peak.
    # We sort the tuple to handle ('a','b') vs ('b','a')
    derived_correct_peak = tuple(sorted(list(unique_major_peaks)[0]))

    # 5. Check the provided answer
    # Find which option corresponds to the derived correct peak
    derived_correct_option = None
    for option_letter, proton_set in options.items():
        if tuple(sorted(list(proton_set))) == derived_correct_peak:
            derived_correct_option = option_letter
            break
    
    if derived_correct_option is None:
        return f"Logic Error: The derived correct peak {derived_correct_peak} does not match any of the options."

    if final_answer_from_llm == derived_correct_option:
        return "Correct"
    else:
        # Provide a detailed reason for the error
        reason = (
            f"The provided answer '{final_answer_from_llm}' is incorrect.\n"
            f"1. Reactants: Maleic anhydride and 1,2,3,4-tetramethyl-1,3-cyclopentadiene.\n"
            f"2. Major Product: Due to severe steric hindrance from the four methyl groups on the diene, the 'exo' adduct is the major product, overriding the standard 'endo' rule.\n"
            f"3. Unique NOESY Peak: In the 'exo' (major) product, the anhydride protons ('anhydride_H', ~3.5 ppm) are spatially close to the vinylic methyl protons ('vinylic_Me', ~1.7 ppm). This interaction is absent in the 'endo' (minor) product.\n"
            f"4. Conclusion: The unique cross-peak connects the signals for 'anhydride_H' and 'vinylic_Me'.\n"
            f"5. Matching Option: This corresponds to option C: A 6H singlet at ~1.7 ppm and a 2H singlet at ~3.5 ppm.\n"
            f"The derived correct option is '{derived_correct_option}', but the provided answer was '{final_answer_from_llm}'."
        )
        return reason

# Run the check
result = check_answer()
print(result)