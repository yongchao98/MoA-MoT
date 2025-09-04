import collections

def check_correctness():
    """
    Checks the correctness of the LLM's answer by analyzing the genetic data.

    The core of the analysis is checking for epistasis:
    - Gene A is epistatic to Gene B if the phenotype of the double mutant (ab)
      is the same as the phenotype of the single mutant of the epistatic gene (a).
    """
    # Experimental data: resistance percentage for each genotype
    phenotypes = {
        "g1": 75,
        "g2": 0,
        "g3": 50,
        "g1g2": 0,
        "g1g3": 10,
        "g2g3": 0,
    }

    # The answer provided by the LLM
    llm_answer = "B"

    # Define a helper function to check for epistasis
    def is_epistatic(epistatic_gene, hypostatic_gene, data):
        """
        Checks if epistatic_gene is epistatic to hypostatic_gene.
        Returns True if the double mutant phenotype matches the epistatic gene's single mutant phenotype.
        """
        # Ensure gene names are in a consistent order for dictionary lookup (e.g., g1g2)
        double_mutant_key = "".join(sorted([epistatic_gene, hypostatic_gene]))
        
        if epistatic_gene not in data or double_mutant_key not in data:
            return None # Data not available to make a conclusion

        return data[double_mutant_key] == data[epistatic_gene]

    # Define the verifiable claims for each option
    # We focus on epistasis as it's directly testable with the given data.
    options_analysis = {
        "A": {
            "claim": "G1 is epistatic towards G3",
            "is_true": is_epistatic("g1", "g3", phenotypes)
        },
        "B": {
            "claim": "G2 is epistatic towards G1",
            "is_true": is_epistatic("g2", "g1", phenotypes)
        },
        "C": {
            "claim": "G1 is epistatic towards G3",
            "is_true": is_epistatic("g1", "g3", phenotypes)
        },
        "D": {
            "claim": "G3 is epistatic towards G1",
            "is_true": is_epistatic("g3", "g1", phenotypes)
        }
    }

    # Evaluate which options are invalidated by the data
    valid_options = []
    rejection_reasons = {}

    for option, analysis in options_analysis.items():
        if analysis["is_true"]:
            valid_options.append(option)
        else:
            # Example: is_epistatic("g1", "g3", phenotypes) -> is 10 == 75? -> False.
            rejection_reasons[option] = f"The claim '{analysis['claim']}' is false based on the data."

    # Final check: Is the LLM's answer the only valid option?
    if llm_answer in valid_options:
        if len(valid_options) == 1:
            return "Correct"
        else:
            return f"Incorrect. The provided answer '{llm_answer}' is consistent with the data, but so are other options: {valid_options}. The question might be ambiguous or the reasoning incomplete."
    else:
        reason = rejection_reasons.get(llm_answer, "The answer is not supported by the data analysis.")
        return f"Incorrect. The provided answer '{llm_answer}' is wrong because a key claim it relies on is false. Reason: {reason}"

# Run the check and print the result
result = check_correctness()
print(result)