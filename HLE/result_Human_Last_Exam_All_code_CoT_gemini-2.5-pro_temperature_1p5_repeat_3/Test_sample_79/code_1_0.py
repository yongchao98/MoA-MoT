import sys
import io

# A helper class to redirect stdout to capture print output for final formatting
class Capturing(list):
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = io.StringIO()
        return self
    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        del self._stringio    # free up some memory
        sys.stdout = self._stdout

def analyze_population_statements():
    """
    Analyzes the statements about a hypothetical population based on a set of given genetic principles.
    """
    # Step 1: Define the conditions given in the problem statement as boolean flags.
    # These represent a population in perfect Hardy-Weinberg equilibrium.
    conditions = {
        "infinite_population": True,  # This implies no genetic drift.
        "random_mating": True,        # Also known as panmixia.
        "no_mutation": True,
        "no_gene_flow": True,
        "no_selection": True,         # Explicitly stated: "no bearing on fitness", "equal fitness".
        "non_overlapping_generations": True
    }

    # Step 2: Evaluate each statement programmatically.
    true_statements = []
    reasoning = []

    # Statement 1: There is no selection occurring on the phenotype measured.
    # Evaluation: This is directly stated in the problem description.
    if conditions["no_selection"]:
        true_statements.append(1)
        reasoning.append(
            "Statement 1 is TRUE.\nThe problem text explicitly states that 'all genotypes have equal fitness' and the trait has 'no bearing on fitness'. This is the definition of no selection."
        )

    # Statement 2: Parents will not raise their offspring.
    # Evaluation: "Non-overlapping generations" means the parental generation is not alive
    # and reproductively active at the same time as the offspring generation. However,
    # this does not strictly forbid all possible forms of parental care (e.g., building a nest,
    # providing care for juveniles before the parents die). Because it is not a logically
    # certain consequence, it does not *always* have to be true.
    reasoning.append(
        "Statement 2 is NOT necessarily true.\n'Non-overlapping generations' is a strong condition, but it does not make all forms of parental care logically impossible. For a statement to 'always be true', there can be no counterexamples. Since a scenario with limited parental care could still fit the description, we cannot confirm this statement as always true."
    )

    # Statement 3: The population will never speciate even in future generations.
    # Evaluation: Speciation is a result of evolution. The fundamental mechanisms of evolution
    # are mutation, selection, genetic drift (ruled out by infinite population), and gene flow.
    # The problem states that all of these are absent. If the conditions persist, the population
    # cannot evolve, and therefore cannot speciate.
    is_evolving = not (conditions["no_mutation"] and \
                       conditions["no_selection"] and \
                       conditions["no_gene_flow"] and \
                       conditions["infinite_population"])
    if not is_evolving:
        true_statements.append(3)
        reasoning.append(
            "Statement 3 is TRUE.\nSpeciation requires evolutionary change. The problem defines a population where all primary mechanisms of evolution (mutation, selection, drift, gene flow) are absent. A non-evolving population cannot speciate."
        )

    # Statement 4: The researcher will not find a substantial difference in the phenotype measured
    # between the west and east groups of individuals.
    # Evaluation: The population is described as a single, randomly mating (panmictic) unit.
    # This means an individual from the west is just as likely to mate with one from the east.
    # This random mixing, along with the absence of localized selection or drift, prevents
    # the formation of any systematic genetic or phenotypic differences between geographic areas.
    # With an infinite sample size, the measurements will be perfectly accurate and show no difference.
    if conditions["random_mating"] and conditions["no_selection"] and conditions["infinite_population"]:
        true_statements.append(4)
        reasoning.append(
            "Statement 4 is TRUE.\nThe population is a single panmictic unit, meaning it mates randomly across the entire region. This prevents the formation of distinct subgroups. Without different selective pressures or genetic drift, there's no mechanism to create a phenotypic difference between the west and east halves."
        )

    # Step 3: Print the analysis and final conclusion.
    print("--- Logical Analysis of Each Statement ---")
    for i, r in enumerate(reasoning):
        print(f"\n[Analysis of Statement {i+1}]")
        print(r)
    print("\n" + "---" * 15)

    print(f"Conclusion: The statements that must always be true are {', '.join(map(str, true_statements))}.")

    # Match the result to the provided answer choices
    options = {
        "A": [1], "B": [2], "C": [3], "D": [4], "E": [1, 2], "F": [1, 3],
        "G": [1, 4], "H": [2, 3], "I": [2, 4], "J": [3, 4], "K": [2, 3, 4],
        "L": [1, 3, 4], "M": [1, 2, 4], "N": [1, 2, 3], "O": [1, 2, 3, 4], "P": []
    }

    final_answer_key = "P" # Default to 'None of the above'
    for key, value in options.items():
        if sorted(value) == sorted(true_statements):
            final_answer_key = key
            break

    print(f"This corresponds to answer choice: {final_answer_key}")
    
    # As requested, output the numbers of the true statements in a final line.
    final_equation = f"Final Tally (True Statements): {true_statements[0]}, {true_statements[1]}, and {true_statements[2]}"
    print(final_equation)
    
    # The final answer in the required format
    print(f"<<<{final_answer_key}>>>")


if __name__ == '__main__':
    analyze_population_statements()