# Define the initial parameters based on the problem description.
self_pollination_rate = 0.05
cross_pollination_rate = 0.95

# Probabilities of genotypes from self-pollination (Wi x Wi)
prob_Wi_from_self = 0.50
prob_ii_from_self = 0.25

# Probabilities of genotypes from cross-pollination (Wi x WW)
prob_Wi_from_cross = 0.50
prob_ii_from_cross = 0.0

# Calculate the contribution of each pollination type to the final offspring population.
# The resistant genotypes are Wi and ii.

# Contribution of Wi from self-pollination
total_freq_Wi_from_self = self_pollination_rate * prob_Wi_from_self
# Contribution of ii from self-pollination
total_freq_ii_from_self = self_pollination_rate * prob_ii_from_self

# Contribution of Wi from cross-pollination
total_freq_Wi_from_cross = cross_pollination_rate * prob_Wi_from_cross
# Contribution of ii from cross-pollination is 0

# Calculate the total frequency of each resistant genotype in the offspring.
total_freq_Wi = total_freq_Wi_from_self + total_freq_Wi_from_cross
total_freq_ii = total_freq_ii_from_self # Since contribution from cross-pollination is 0

# The total frequency of resistant offspring is the sum of the frequencies of Wi and ii genotypes.
total_resistant_freq = total_freq_Wi + total_freq_ii

# Convert the final frequency to a percentage.
total_resistant_percentage = total_resistant_freq * 100

# Print the step-by-step calculation.
print(f"The parent plant (Wi) undergoes 5% self-pollination and 95% cross-pollination with the WW population.")
print(f"Resistant phenotypes are conferred by Wi and ii genotypes.\n")
print(f"1. Calculate the total frequency of the resistant Wi genotype:")
print(f"   - From self-pollination: {self_pollination_rate} (rate) * {prob_Wi_from_self} (prob) = {total_freq_Wi_from_self}")
print(f"   - From cross-pollination: {cross_pollination_rate} (rate) * {prob_Wi_from_cross} (prob) = {total_freq_Wi_from_cross}")
print(f"   - Total Freq(Wi) = {total_freq_Wi_from_self} + {total_freq_Wi_from_cross} = {total_freq_Wi}\n")

print(f"2. Calculate the total frequency of the resistant ii genotype:")
print(f"   - From self-pollination: {self_pollination_rate} (rate) * {prob_ii_from_self} (prob) = {total_freq_ii_from_self}")
print(f"   - From cross-pollination: The ii genotype cannot be produced.")
print(f"   - Total Freq(ii) = {total_freq_ii_from_self}\n")

print(f"3. Calculate the total frequency of resistant offspring:")
print(f"   - Total Resistant Freq = Freq(Wi) + Freq(ii)")
print(f"   - Total Resistant Freq = {total_freq_Wi} + {total_freq_ii} = {total_resistant_freq}\n")

print(f"Final Answer: The theoretical percentage of drought-resistant offspring is {total_resistant_freq} * 100 = {total_resistant_percentage}%.")