def analyze_allele_detection(allele_frequency, sequencing_error_rate, depth):
    """
    Calculates the expected number of true allele reads vs. error reads.
    """
    expected_true_reads = depth * allele_frequency
    # We multiply by (1 - allele_frequency) as errors can't happen on a true allele read
    expected_error_reads = depth * (1 - allele_frequency) * sequencing_error_rate
    return expected_true_reads, expected_error_reads

def main():
    """
    Main function to run the analysis and print results.
    """
    # Define parameters for the simulation
    low_freq_allele = 0.01  # A hypothetical low-frequency allele at 1%
    error_rate = 0.001        # A typical high-quality sequencing error rate of 0.1%
    low_depth = 40
    high_depth = 80

    print("This script demonstrates why increasing sequencing depth is critical for detecting low-frequency alleles.")
    print(f"We will model a scenario with a true allele at {low_freq_allele*100}% frequency and a sequencing error rate of {error_rate*100}%.")
    print("-" * 50)

    # --- Calculation for Low Depth ---
    print(f"Analysis at {low_depth}X Depth:")
    true_reads_low, error_reads_low = analyze_allele_detection(low_freq_allele, error_rate, low_depth)
    print("Expected number of TRUE allele reads:")
    print(f"{low_depth} (depth) * {low_freq_allele} (frequency) = {true_reads_low:.2f}")
    print("Expected number of reads from sequencing ERROR:")
    print(f"{low_depth} (depth) * (1 - {low_freq_allele}) * {error_rate} (error rate) = {error_reads_low:.3f}")
    print("-" * 50)

    # --- Calculation for High Depth ---
    print(f"Analysis at {high_depth}X Depth (as proposed in Option A):")
    true_reads_high, error_reads_high = analyze_allele_detection(low_freq_allele, error_rate, high_depth)
    print("Expected number of TRUE allele reads:")
    print(f"{high_depth} (depth) * {low_freq_allele} (frequency) = {true_reads_high:.2f}")
    print("Expected number of reads from sequencing ERROR:")
    print(f"{high_depth} (depth) * (1 - {low_freq_allele}) * {error_rate} (error rate) = {error_reads_high:.3f}")
    print("-" * 50)
    
    print("Conclusion: Increasing depth from 40X to 80X doubles the expected signal (from "
          f"{true_reads_low:.2f} to {true_reads_high:.2f} reads),")
    print("providing much stronger statistical evidence to distinguish a real allele from background noise.")
    print("This, combined with a 'thorough library preparation' to minimize error and bias, makes Option A the best choice.")

if __name__ == "__main__":
    main()