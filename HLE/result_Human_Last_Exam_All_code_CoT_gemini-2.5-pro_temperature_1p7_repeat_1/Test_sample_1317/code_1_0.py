import sys
import io

# Backup original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

def solve():
    """
    This function demonstrates the relationship between gene transfer time and
    recombination frequency in an E. coli interrupted mating experiment.
    """
    # Step 1: Define the gene order and their relative times of entry.
    # The problem states the gene order is thr-azi-gal and thr+ is transferred first.
    # This means the origin of transfer (oriT) is located just before the thr gene.
    # We assign arbitrary 'time of entry' units, where smaller numbers mean earlier entry.
    gene_entry_times = {
        'thr': 2,   # Enters first (earliest time)
        'azy': 9,   # Enters second
        'gal': 17   # Enters third (latest time)
    }

    print("--- Interrupted Mating Experiment Simulation ---")
    print("\nGiven Gene Order on Chromosome: oriT -> thr -> azi -> gal")
    print("Assigned Relative Time of Entry (e.g., in minutes):")
    for gene, time in gene_entry_times.items():
        print(f"- Gene '{gene}': {time} minutes")

    # Step 2: Model the recombination frequency.
    # The frequency of recombination is inversely proportional to the time of entry.
    # Genes that enter earlier (smaller time) have a higher chance of being
    # transferred before the mating is interrupted. We use a simple inverse formula
    # to model this relationship.
    print("\nCalculating Relative Recombination Frequencies:")
    print("Formula: Frequency = 100 / Time of Entry")

    recombination_frequencies = {}
    for gene, time in gene_entry_times.items():
        # Using 100 as a base for easier visualization
        frequency = 100 / time
        recombination_frequencies[gene] = frequency
        # Output the calculation for clarity
        print(f"Frequency for '{gene}': 100 / {time} = {frequency:.2f}")

    # Step 3: Identify the gene with the highest frequency and its location.
    highest_freq_gene = max(recombination_frequencies, key=recombination_frequencies.get)
    highest_freq_value = recombination_frequencies[highest_freq_gene]

    print(f"\nResult: The highest recombination frequency ({highest_freq_value:.2f}) is observed for the '{highest_freq_gene}' gene.")

    # Step 4: Draw the final conclusion.
    print("\nConclusion:")
    print(f"The highest frequency of recombinants occurs for the gene that is transferred first ('{highest_freq_gene}').")
    print("This gene is located closest to the origin of transfer (oriT).")
    print("Therefore, the genetic location where the gradient of recombination is highest is the origin itself.")
    print("This corresponds to the position 'Immediately before thr+'.")

solve()

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()
print(output)
print("<<<D>>>")