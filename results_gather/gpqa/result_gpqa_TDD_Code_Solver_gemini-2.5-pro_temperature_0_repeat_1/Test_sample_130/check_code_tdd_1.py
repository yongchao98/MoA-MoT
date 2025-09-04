import re

def check_diels_alder_noesy_answer():
    """
    This function checks the correctness of the provided answer to a chemistry problem
    involving a Diels-Alder reaction and NOESY NMR analysis.

    It codifies the chemical principles and logical steps to derive the correct answer
    and compares it with the provided one.
    """
    try:
        # --- Part 1: Define Chemical Principles and System ---

        # Principle 1: In a Diels-Alder reaction, the 'endo' product is typically the major,
        # kinetically-favored product. The 'exo' is the minor product.
        major_isomer = "endo"
        minor_isomer = "exo"

        # Principle 2: Define the key proton groups in the product and their expected signals
        # based on the question's options and general chemical knowledge.
        proton_groups = {
            "anhydride_H":   {"protons": 2, "shift": 3.5},
            "vinylic_Me":    {"protons": 6, "shift": 1.7},
            "bridgehead_Me": {"protons": 6, "shift": 1.0},
            "bridge_H":      {"protons": 1, "shift": 1.5} # Represents one of the two bridge protons
        }

        # Principle 3: Define spatial proximities (<5 Å) for NOESY analysis in each isomer.
        # This is the core chemical insight.
        # Key: a sorted tuple of proton group names. Value: 'close' or 'far'.
        proximities = {
            "endo": { # Major product
                tuple(sorted(("anhydride_H", "vinylic_Me"))): "close",
                tuple(sorted(("anhydride_H", "bridgehead_Me"))): "far",
            },
            "exo": { # Minor product
                tuple(sorted(("anhydride_H", "vinylic_Me"))): "far",
                tuple(sorted(("anhydride_H", "bridgehead_Me"))): "close",
            }
        }

        # --- Part 2: Logical Deduction ---

        # The question asks for a cross-peak present in the major product but absent in the minor.
        # This means we need a pair of protons that is 'close' in the 'endo' isomer
        # and 'far' in the 'exo' isomer.
        unique_major_correlation = None
        for pair, distance_in_major in proximities[major_isomer].items():
            if distance_in_major == "close" and proximities[minor_isomer].get(pair) == "far":
                unique_major_correlation = pair
                break

        if not unique_major_correlation:
            return "Logic Error: The chemical model did not yield a unique correlation for the major product."

        # --- Part 3: Match to Options ---

        options = {
            "A": ["A 6H singlet at ~1 ppm", "a 1H doublet at ~1.5 ppm"],
            "B": ["A 6H singlet at ~1.7 ppm", "a 2H singlet at ~3.5 ppm"],
            "C": ["A 1H doublet at ~1.5 ppm", "a 2H singlet at ~3.5 ppm"],
            "D": ["A 6H singlet at ~1 ppm", "a 6H singlet at ~1.7 ppm"]
        }

        def get_group_from_desc(desc_text):
            match = re.search(r"(\d+)H.*?([\d\.]+)\s*ppm", desc_text)
            if not match: return None
            protons, shift = int(match.group(1)), float(match.group(2))
            for name, details in proton_groups.items():
                if abs(details["shift"] - shift) < 0.2 and details["protons"] == protons:
                    return name
            return None

        derived_correct_option = None
        for letter, descriptions in options.items():
            group1_name = get_group_from_desc(descriptions[0])
            group2_name = get_group_from_desc(descriptions[1])
            if group1_name and group2_name:
                if tuple(sorted((group1_name, group2_name))) == unique_major_correlation:
                    derived_correct_option = letter
                    break

        if not derived_correct_option:
            return f"Logic Error: The correctly identified proton pair {unique_major_correlation} could not be matched to any of the options."

        # --- Part 4: Final Verification ---
        llm_answer = "B"
        llm_reasoning = """
        An analysis of the reaction and the resulting stereoisomers is required to solve this problem.
        1.  Reactants and Reaction: The reaction is a Diels-Alder cycloaddition between maleic anhydride (the dienophile) and 1,2,3,4-tetramethyl-1,3-cyclopentadiene (the diene).
        2.  Stereochemistry (Products): ... Endo Product (Major): ... anhydride ring is tucked "under" ... Exo Product (Minor): ... anhydride ring is oriented "away" ...
        3.  NOESY Spectroscopy: ... detects protons that are close to each other in space ...
        4.  Structural Analysis and Proton Assignment: ... Anhydride Protons (2H): ... ~3.5 ppm. ... Vinylic Methyl Protons (6H): ... ~1.7 ppm. ...
        5.  Proximity in Isomers: ... Major (Endo) Product: ... anhydride protons (~3.5 ppm) very close to the vinylic methyl groups (~1.7 ppm). ... Minor (Exo) Product: ... anhydride protons (~3.5 ppm) are far from the vinylic methyl groups (~1.7 ppm) ...
        6.  Conclusion: The NOE cross-peak ... is the one connecting the anhydride protons and the vinylic methyl protons. ... This corresponds to option B.
        """

        if llm_answer == derived_correct_option:
            # Check for key logical steps in the reasoning
            if "Endo Product (Major)" in llm_reasoning and \
               "anhydride protons (~3.5 ppm) very close to the vinylic methyl groups (~1.7 ppm)" in llm_reasoning and \
               "anhydride protons (~3.5 ppm) are far from the vinylic methyl groups (~1.7 ppm)" in llm_reasoning:
                return "Correct"
            else:
                return "Incorrect: The final answer letter is correct, but the provided reasoning is flawed or incomplete."
        else:
            return f"Incorrect: The provided answer is {llm_answer}, but the correct answer based on stereochemical analysis is {derived_correct_option}. The unique NOESY cross-peak in the major (endo) product is between the anhydride protons (~3.5 ppm) and the vinylic methyl protons (~1.7 ppm), as these are brought into close proximity in the endo isomer but are far apart in the exo isomer."

    except Exception as e:
        return f"An error occurred during the checking process: {e}"

# Run the check
result = check_diels_alder_noesy_answer()
print(result)