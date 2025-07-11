import random

def run_simulation():
    # --- Simulation Parameters ---
    POP_SIZE = 100
    GENERATIONS = 150
    NUM_GENES = 10
    MUTATION_RATE = 0.02  # Chance per generation a mutation hits a specific gene/gene family

    # --- Scenario 1: Single-Copy Genes ---
    # Genome is a list of genes. 1 = functional, 0 = non-functional.
    # e.g., [1, 1, 1, 1, 1]
    pop_single_copy = [[1] * NUM_GENES for _ in range(POP_SIZE)]

    for _ in range(GENERATIONS):
        # Asexual reproduction: each new individual is a copy of a parent
        for i in range(POP_SIZE):
            # A mutation may hit one gene, breaking it permanently.
            if random.random() < MUTATION_RATE * NUM_GENES:
                gene_to_mutate = random.randint(0, NUM_GENES - 1)
                pop_single_copy[i][gene_to_mutate] = 0

    # --- Scenario 2: Multigene Families ---
    # Genome has gene families with 3 copies each.
    # e.g., [[1,1,1], [1,1,1], [1,1,1]]
    num_copies_per_family = 3
    pop_multi_gene = [[[1] * num_copies_per_family for _ in range(NUM_GENES)] for _ in range(POP_SIZE)]

    for _ in range(GENERATIONS):
        for i in range(POP_SIZE):
            # A mutation may hit one gene family.
            if random.random() < MUTATION_RATE * NUM_GENES:
                family_to_mutate = random.randint(0, NUM_GENES - 1)
                # It breaks one copy within the family.
                copy_to_mutate = random.randint(0, num_copies_per_family - 1)
                pop_multi_gene[i][family_to_mutate][copy_to_mutate] = 0

    # --- Analysis & Output ---
    # Calculate how many individuals have lost at least one essential gene function.

    # For single-copy genomes, loss of function happens if any gene is 0.
    lost_function_single = 0
    for genome in pop_single_copy:
        if 0 in genome:
            lost_function_single += 1

    # For multigene families, loss of function only happens if ALL copies in a family are 0.
    lost_function_multi = 0
    for genome in pop_multi_gene:
        for gene_family in genome:
            if sum(gene_family) == 0:
                lost_function_multi += 1
                break  # This individual has lost a function, count it and move to the next.

    print("SIMULATION: Testing compensatory mechanisms for genetic deterioration.")
    print("-" * 60)
    print(f"In a population with limited recombination, deleterious mutations accumulate (Muller's Ratchet).")
    print(f"We simulated this for {GENERATIONS} generations in a population of {POP_SIZE}.\n")

    print("--- Scenario A: Genomes with Single-Copy Genes ---")
    print("A single mutation in a gene leads to permanent loss of its function.")
    print(f"Result: {lost_function_single} out of {POP_SIZE} individuals suffered a complete loss of at least one gene function.\n")

    print("--- Scenario D: Genomes with Multigene Families ---")
    print("A gene's function is only lost if all copies in its family are mutated.")
    print(f"Result: {lost_function_multi} out of {POP_SIZE} individuals suffered a complete loss of at least one gene function.\n")

    print("-" * 60)
    print("CONCLUSION:")
    print("The simulation shows that having redundant gene copies (Multigene Families) provides a powerful buffer against genetic deterioration.")
    print("This strongly supports Multigene Families as the intrinsic architectural feature that compensates for limited recombination.")
    print("The correct choice is D.")

run_simulation()
<<<D>>>