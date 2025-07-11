import random

def calculate_selection_coefficient(test_fitness, reference_fitness):
    """
    Calculates the selection coefficient (s).
    s = 1 - (w_test / w_ref)
    A positive 's' indicates the test strain is less fit than the reference.
    A negative 's' indicates the test strain is more fit (e.g., hybrid vigor).
    """
    if reference_fitness == 0:
        return float('inf')
    return 1 - (test_fitness / reference_fitness)

def run_experiment_simulation():
    """
    Simulates an experiment to measure the cost of gene flow in yeast.
    This follows the logic of comparing no-gene-flow lines to hybrids
    and accounting for the effects of meiosis.
    """
    # Step 1: Define the "no gene flow" parental lines.
    # Assume they are well-adapted and have a standard fitness of 1.0.
    parent_1_fitness = 1.0
    parent_2_fitness = 1.0
    no_gene_flow_reference_fitness = (parent_1_fitness + parent_2_fitness) / 2

    print("--- Experimental Setup ---")
    print(f"Fitness of Parent 1 (No Gene Flow Line): {parent_1_fitness}")
    print(f"Fitness of Parent 2 (No Gene Flow Line): {parent_2_fitness}")
    print(f"Reference Fitness (Average of Parents): {no_gene_flow_reference_fitness}\n")

    # Step 2: Measure the fitness of F1 hybrids.
    # Gene flow creates F1 hybrids. They may show "hybrid vigor" (heterosis).
    f1_hybrid_fitness = 1.05  # Assume a 5% fitness advantage for F1
    s_f1 = calculate_selection_coefficient(f1_hybrid_fitness, no_gene_flow_reference_fitness)

    print("--- F1 Hybrid Generation (Initial Gene Flow) ---")
    print(f"Fitness of F1 Hybrid: {f1_hybrid_fitness}")
    print("s_f1 = 1 - (F1_Fitness / Reference_Fitness)")
    print(f"s_f1 = 1 - ({f1_hybrid_fitness} / {no_gene_flow_reference_fitness}) = {s_f1:.3f}")
    print("The negative selection coefficient suggests F1 hybrids are more fit (hybrid vigor).\n")

    # Step 3: Check for effects of meiosis by creating an F2 generation.
    # F1 hybrids undergo meiosis and mate ("within mating").
    # Recombination can break up co-adapted gene complexes, causing outbreeding depression.
    # Let's simulate a drop in fitness in the F2 generation.
    f2_generation_fitness = 0.80 # Assume a 20% fitness cost revealed after meiosis
    s_f2 = calculate_selection_coefficient(f2_generation_fitness, no_gene_flow_reference_fitness)

    print("--- F2 Generation (After Meiosis and Within Mating) ---")
    print(f"Average Fitness of F2 Generation: {f2_generation_fitness}")
    print("s_f2 = 1 - (F2_Fitness / Reference_Fitness)")
    print(f"s_f2 = 1 - ({f2_generation_fitness} / {no_gene_flow_reference_fitness}) = {s_f2:.3f}")
    print("The positive selection coefficient reveals a significant fitness cost (outbreeding depression) after meiosis.\n")

    print("--- Conclusion ---")
    print("A complete measurement of the cost of gene flow requires calculating the selection coefficient for both F1 and F2 generations.")
    print("This allows for the quantification of costs hidden in the F1 generation that are only revealed by the effects of meiosis.")

# Run the simulation
run_experiment_simulation()