import collections

def check_diels_alder_noesy():
    """
    Checks the correctness of the answer for the Diels-Alder NOESY problem.

    The logic follows these steps:
    1.  Define the proton types and their expected NMR signals based on the problem description.
    2.  Define the answer options in terms of these proton types.
    3.  Establish the stereochemical rules for the reaction:
        a. A Diels-Alder reaction produces 'endo' and 'exo' products.
        b. Crucially, due to severe steric hindrance from the tetramethylated diene, the 'exo' product is major, reversing the usual Alder-endo rule.
    4.  Define the spatial proximities (< 5 Å) for both 'endo' and 'exo' isomers, which determine NOESY cross-peaks.
    5.  Identify the NOESY cross-peak that is unique to the major product ('exo') but absent in the minor product ('endo').
    6.  Determine which answer option corresponds to this unique cross-peak.
    7.  Compare the derived correct option with the provided answer ('C').
    """
    
    # Step 1: Define proton types and their assigned NMR signals
    protons = {
        "H_anhydride": "2H singlet at ~3.5 ppm",
        "Me_vinyl": "6H singlet at ~1.7 ppm",
        "H_bridge": "1H doublet at ~1.5 ppm",
        "Me_bridgehead": "6H singlet at ~1.0 ppm"
    }

    # Step 2: Define the options based on the proton types
    # Using frozenset to make the pair order-independent
    options = {
        "A": frozenset(["H_bridge", "H_anhydride"]),
        "B": frozenset(["Me_bridgehead", "H_bridge"]),
        "C": frozenset(["Me_vinyl", "H_anhydride"]),
        "D": frozenset(["Me_bridgehead", "Me_vinyl"])
    }
    
    provided_answer = "C"

    # Step 3: Establish stereochemical rules
    # The key insight is that severe steric hindrance from the four methyl groups
    # on the diene reverses the usual Alder-endo rule.
    major_product = "exo"
    minor_product = "endo"
    
    reasoning_log = []
    reasoning_log.append(f"Step 1: The reaction is a Diels-Alder between maleic anhydride and 1,2,3,4-tetramethyl-1,3-cyclopentadiene.")
    reasoning_log.append(f"Step 2: The reaction produces two stereoisomers: '{major_product}' and '{minor_product}'.")
    reasoning_log.append(f"Step 3: Due to severe steric hindrance, the '{major_product}' adduct is the major product, reversing the typical Alder-endo rule.")

    # Step 4: Define spatial proximities for NOESY cross-peaks in each isomer
    # In the 'exo' adduct, the anhydride protons are on the endo face, close to the vinylic methyls.
    exo_proximities = {
        frozenset(["H_anhydride", "Me_vinyl"]) 
        # also close to syn-H_bridge, but this is the most distinguishing interaction
    }
    
    # In the 'endo' adduct, the anhydride protons are on the exo face, close to the bridgehead methyls and anti-H_bridge.
    endo_proximities = {
        frozenset(["H_anhydride", "Me_bridgehead"]),
        frozenset(["H_anhydride", "H_bridge"]) # This represents the anti-H_bridge
    }
    
    reasoning_log.append(f"Step 4: Analyzing spatial proximities for NOESY:")
    reasoning_log.append(f" - In the major ('{major_product}') product, H_anhydride is close to Me_vinyl.")
    reasoning_log.append(f" - In the minor ('{minor_product}') product, H_anhydride is far from Me_vinyl.")

    # Step 5: Find the cross-peak present in the major product but absent in the minor
    unique_major_cross_peaks = exo_proximities - endo_proximities
    
    if not unique_major_cross_peaks:
        return "Logic Error: Could not find a unique cross-peak for the major product based on the defined proximities."
        
    # We expect one key distinguishing cross-peak
    distinguishing_peak_protons = list(unique_major_cross_peaks)[0]
    
    reasoning_log.append(f"Step 5: The cross-peak present in the major product but absent in the minor must connect the protons: {list(distinguishing_peak_protons)}.")

    # Step 6: Find which option corresponds to this cross-peak
    correct_option = None
    for option, proton_pair in options.items():
        if proton_pair == distinguishing_peak_protons:
            correct_option = option
            break
            
    if correct_option is None:
        return "Logic Error: The derived correct proton pair does not match any of the options."

    p1, p2 = list(distinguishing_peak_protons)
    reasoning_log.append(f"Step 6: This corresponds to the interaction between the '{protons[p1]}' and the '{protons[p2]}'. This is option {correct_option}.")

    # Step 7: Compare with the provided answer
    if provided_answer == correct_option:
        return "Correct"
    else:
        reasoning_log.append(f"\nError: The provided answer is '{provided_answer}', but the correct answer derived from chemical principles is '{correct_option}'.")
        reasoning_log.append(f"Reason: The provided answer does not correspond to the unique NOESY interaction expected for the major ('{major_product}') product.")
        return "\n".join(reasoning_log)

# Run the check
result = check_diels_alder_noesy()
print(result)