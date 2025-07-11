import random

def simulate_gene_duplication_fates():
    """
    This script conceptually models the primary fates of a duplicated gene
    to illustrate why subfunctionalization and neofunctionalization are
    key mechanisms for retention and divergence.
    """

    # 1. Define an ancestral gene with multiple functions (subfunctions)
    ancestral_gene = {
        "name": "AncestralGene-XYZ",
        "functions": ["catalyzes reaction A", "binds to protein B"]
    }
    print("--- Initial State ---")
    print(f"Original Gene: {ancestral_gene['name']}")
    print(f"Functions: {', '.join(ancestral_gene['functions'])}\n")

    print("--- Gene Duplication Occurs ---\n")

    # FATE 1: PSEUDOGENIZATION (The most common fate)
    print("--- Scenario 1: Pseudogenization ---")
    copy1 = ancestral_gene.copy()
    copy2 = {"name": "AncestralGene-XYZ-copy", "functions": []} # Loses function
    print(f"Gene 1 ({copy1['name']}) retains its function.")
    print(f"Gene 2 ({copy2['name']}) accumulates mutations and becomes a non-functional pseudogene.")
    print("Result: Only one functional gene copy is retained. No functional divergence.\n")


    # FATE 2: NEOFUNCTIONALIZATION (Acquiring a new function)
    print("--- Scenario 2: Neofunctionalization ---")
    copy1 = ancestral_gene.copy()
    copy2 = {
        "name": "AncestralGene-XYZ-copy",
        "functions": ["regulates pathway C"] # A completely new function
    }
    print(f"Gene 1 ({copy1['name']}) retains original functions: {', '.join(copy1['functions'])}")
    print(f"Gene 2 ({copy2['name']}) evolves a new function: {', '.join(copy2['functions'])}")
    print("Result: Both gene copies are retained and have diverged. This creates evolutionary novelty.\n")


    # FATE 3: SUBFUNCTIONALIZATION (Partitioning of functions)
    print("--- Scenario 3: Subfunctionalization ---")
    copy1 = {
        "name": "AncestralGene-XYZ-alpha",
        "functions": [ancestral_gene['functions'][0]] # Keeps only the first function
    }
    copy2 = {
        "name": "AncestralGene-XYZ-beta",
        "functions": [ancestral_gene['functions'][1]] # Keeps only the second function
    }
    print("The ancestral functions are partitioned between the two copies:")
    print(f"Gene 1 ({copy1['name']}) now specializes in: {', '.join(copy1['functions'])}")
    print(f"Gene 2 ({copy2['name']}) now specializes in: {', '.join(copy2['functions'])}")
    print("Result: Both copies must be retained to perform the full ancestral role. They have diverged through specialization.")

simulate_gene_duplication_fates()