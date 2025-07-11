import numpy as np

def simulate_allele_detection():
    """
    This script simulates the detection of a low-frequency allele in a sequencing experiment.
    It demonstrates how increasing sequencing depth improves the probability of
    confidently detecting the allele above the background noise of sequencing errors.
    """

    # --- Simulation Parameters ---
    # The actual frequency of the rare drug-resistance allele in the microbial population.
    # Let's say it's very low, 0.5% or 0.005.
    true_allele_frequency = 0.005

    # Per-base sequencing error rate. For high-quality Illumina reads, this is very low, e.g., 0.1%.
    # This is the rate at which a reference base is misread as our specific allele of interest.
    sequencing_error_rate = 0.001

    # The minimum number of reads required to support the allele before we believe it's real
    # and not just a random sequencing error. A common threshold is 3 or more reads.
    detection_threshold = 3

    # The sequencing depths we want to compare, as mentioned in the question.
    # We add higher depths to show the effect more clearly.
    depths_to_test = [40, 80, 500, 1000]

    # Number of times we run the simulation for each depth to get a stable probability.
    num_simulations = 10000

    print("--- Simulating Detection of a Low-Frequency Allele ---")
    print(f"True Allele Frequency: {true_allele_frequency * 100}%")
    print(f"Sequencing Error Rate: {sequencing_error_rate * 100}%")
    print(f"Detection requires at least {detection_threshold} reads with the allele.")
    print("-" * 50)

    # --- Run Simulation for each depth ---
    for depth in depths_to_test:
        detection_count = 0
        for _ in range(num_simulations):
            # Step 1: From the total 'depth' reads, how many are from cells that truly have the rare allele?
            # We use a binomial distribution to model this sampling process.
            true_alt_reads = np.random.binomial(n=depth, p=true_allele_frequency)
            true_ref_reads = depth - true_alt_reads

            # Step 2: Account for sequencing errors that can mimic the allele.
            # How many of the true 'ref' reads are incorrectly sequenced as the 'alt' allele (False Positives)?
            false_positives = np.random.binomial(n=true_ref_reads, p=sequencing_error_rate)

            # The total number of OBSERVED reads with the allele is the sum of true reads and error-generated reads.
            # (For simplicity, we assume all true_alt_reads are sequenced correctly).
            observed_alt_reads = true_alt_reads + false_positives

            # Step 3: Check if the observed count meets our detection threshold.
            if observed_alt_reads >= detection_threshold:
                detection_count += 1

        # Calculate the overall probability of detection at this depth
        detection_probability = detection_count / num_simulations

        print(f"At {depth}X sequencing depth:")
        print(f"  Probability of detecting the allele: {detection_probability * 100:.2f}%")

    print("-" * 50)
    print("Simulation Conclusion: Increasing sequencing depth dramatically increases the")
    print("probability of detecting a low-frequency allele above the noise threshold.")
    print("This demonstrates why increasing depth from 40X to 80X, as suggested in")
    print("the problem, is a critical step in the correct strategy.")

# Run the simulation
simulate_allele_detection()