import math

# The problem describes an interrupted mating experiment. In this process, genes are
# transferred linearly from a donor to a recipient cell. The key principle is that the
# connection between the cells is fragile and can break at any time. This means that
# genes transferred earlier have a much higher chance of making it into the recipient
# than genes transferred later. The frequency of recombinants for a specific gene is
# directly related to this transfer probability.

# --- Step 1: Define the Model Parameters ---
# The gene transfer order is given as thr -> azi -> gal.
# We can assign arbitrary entry times to represent this sequence. The absolute
# values are not important, only their relative order matters for the outcome.
gene_entry_times = {
    'thr': 10,  # First gene, enters at 10 minutes
    'azi': 20,  # Second gene, enters at 20 minutes
    'gal': 30   # Third gene, enters at 30 minutes
}

# We need a value to represent the fragility of the mating connection.
# Let's assume a constant probability that the connection breaks per unit of time.
prob_interruption_per_unit_time = 0.05
# From this, we can find the probability that the connection does NOT break.
prob_no_interruption_per_unit_time = 1 - prob_interruption_per_unit_time

print("Principle: The probability of a gene being transferred successfully decreases as its time of entry increases.")
print("This can be modeled with the equation:")
print("  P(transfer) = (P(no interruption per unit time)) ^ (time of entry)\n")

highest_frequency = -1.0
location_of_highest_frequency = None

# --- Step 2: Calculate Transfer Probability for Each Gene ---
print("--- Calculating Relative Transfer Frequencies ---")
for gene, time in sorted(gene_entry_times.items(), key=lambda item: item[1]):
    # The probability of successful transfer is the chance that the connection
    # survives for the entire duration of the gene's transfer time.
    transfer_probability = prob_no_interruption_per_unit_time ** time
    
    # Print the equation with the specific numbers for this gene.
    print(f"For gene '{gene}' (enters at time={time}):")
    print(f"  P(transfer) = {prob_no_interruption_per_unit_time:.2f} ^ {time} = {transfer_probability:.4f}")
    
    # Keep track of the highest frequency found so far.
    if transfer_probability > highest_frequency:
        highest_frequency = transfer_probability
        location_of_highest_frequency = gene

# --- Step 3: State the Conclusion ---
print("\n--- Conclusion ---")
print(f"The simulation shows the highest transfer frequency ({highest_frequency:.4f}) is for the '{location_of_highest_frequency}' gene.")
print("This is because it is the first genetic marker transferred from the donor to the recipient.")
print("Therefore, the highest frequency of recombinants will be observed at the genetic location of the first transferred gene, which is 'thr'.")
print("Among the given choices, the location 'Immediately after thr+' best describes the position of this first marker relative to the start of the transfer process.")
