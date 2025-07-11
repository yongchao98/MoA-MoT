import textwrap

def analyze_population_scenario():
    """
    Analyzes the population genetics scenario to determine which statements must be true.
    """
    print("Step-by-step analysis of each statement:\n")

    # --- Statement 1 ---
    s1_analysis = """
    Statement 1: There is no selection occurring on the phenotype measured.

    Analysis: The problem explicitly states that 'all genotypes have equal fitness' and the phenotype 'has no bearing on fitness'. This is the definition of no natural selection.
    
    Conclusion: Statement 1 MUST be true.
    """
    print(textwrap.dedent(s1_analysis))

    # --- Statement 2 ---
    s2_analysis = """
    Statement 2: Parents will not raise their offspring.

    Analysis: The problem states that generations are 'discrete and non-overlapping'. This is a modeling assumption meaning that the parental generation dies before the offspring generation reaches reproductive age. It prevents mating between generations. However, it does not rule out parental care. For example, parents could raise their young to independence and then die. Since parental care is not explicitly ruled out, we cannot say this statement must always be true.

    Conclusion: Statement 2 does not have to be true.
    """
    print(textwrap.dedent(s2_analysis))

    # --- Statement 3 ---
    s3_analysis = """
    Statement 3: The population will never speciate even in future generations.

    Analysis: The current conditions (no mutation, no selection, no drift, no gene flow, random mating) describe a population in Hardy-Weinberg Equilibrium, which is not evolving. However, these conditions are only described as 'currently true'. They are not guaranteed to hold for all future generations. A future mutation, a new selective pressure, or the formation of a geographic barrier could disrupt the equilibrium and lead to speciation. The word 'never' is too strong an assertion.

    Conclusion: Statement 3 does not have to be true.
    """
    print(textwrap.dedent(s3_analysis))

    # --- Statement 4 ---
    s4_analysis = """
    Statement 4: The researcher will not find a substantial difference in the phenotype measured between the west and east groups of individuals.

    Analysis: The population is described as 'randomly mating'. This implies that mating occurs without regard to geographic location within the region. Such panmixia ensures that allele frequencies are homogenized throughout the entire population. Therefore, the east and west halves are simply arbitrary samples from a single, genetically uniform population. Since the population is infinite and the measurements are perfect, the measured phenotypic distributions of the two halves will be identical.

    Conclusion: Statement 4 MUST be true.
    """
    print(textwrap.dedent(s4_analysis))

    # --- Final Conclusion ---
    final_conclusion = """
    Summary:
    - Statement 1 must be true.
    - Statement 2 is not necessarily true.
    - Statement 3 is not necessarily true.
    - Statement 4 must be true.

    Therefore, the combination of statements that must always be true is 1 and 4.
    This corresponds to answer choice G.
    """
    print("-----------------------------------------")
    print(textwrap.dedent(final_conclusion))


if __name__ == '__main__':
    analyze_population_scenario()
    # The final answer is determined by the analysis above.
    # Statements 1 and 4 must be true.
    final_answer = "G"
    print(f"<<<{final_answer}>>>")