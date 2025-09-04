def check_genetic_conclusion():
    """
    Checks the correctness of the provided answer by analyzing the experimental data
    based on principles of genetic interaction.
    """
    # Resistance percentages for each mutant genotype
    phenotypes = {
        "wt": 100,
        "g1": 75,
        "g2": 0,
        "g3": 50,
        "g1g2": 0,
        "g1g3": 10,
        "g2g3": 0,
    }

    # The provided answer is C
    # Statement C: "G2 is a transcription factor, G1 and G3 show gene redundancy, G1 is epistatic towards G3"

    # Part 1: Check if G2 is a transcription factor (i.e., epistatic to others)
    is_g2_epistatic_to_g1 = (phenotypes["g1g2"] == phenotypes["g2"])
    is_g2_epistatic_to_g3 = (phenotypes["g2g3"] == phenotypes["g2"])
    is_g2_transcription_factor = is_g2_epistatic_to_g1 and is_g2_epistatic_to_g3

    # Part 2: Check if G1 and G3 show gene redundancy (synergistic interaction)
    # A synergistic interaction means the double mutant phenotype is more severe (lower resistance)
    # than either single mutant phenotype.
    is_synergistic = (phenotypes["g1g3"] < phenotypes["g1"]) and (phenotypes["g1g3"] < phenotypes["g3"])
    is_redundant = is_synergistic

    # Part 3: Check if G1 is epistatic towards G3
    # This would mean the g1g3 phenotype is the same as the g1 phenotype.
    is_g1_epistatic_to_g3 = (phenotypes["g1g3"] == phenotypes["g1"])

    # Evaluate the full statement from answer C
    if not is_g2_transcription_factor:
        return "Incorrect. The claim that 'G2 is a transcription factor' is not supported. G2's knockout phenotype (0%) is not epistatic to both G1 and G3 based on the double mutant data."
    
    if not is_redundant:
        return "Incorrect. The claim that 'G1 and G3 show gene redundancy' is not supported. The g1g3 double mutant phenotype (10%) is not more severe than both the g1 (75%) and g3 (50%) single mutant phenotypes."

    if not is_g1_epistatic_to_g3:
        return f"Incorrect. The answer C claims that 'G1 is epistatic towards G3'. This is false. For G1 to be epistatic to G3, the phenotype of the g1g3 double mutant ({phenotypes['g1g3']}%) should be the same as the g1 single mutant ({phenotypes['g1']}%). The data clearly shows they are different."

    # If all parts of the statement were true, the answer would be correct.
    return "Correct"

# Execute the check and print the result.
result = check_genetic_conclusion()
print(result)