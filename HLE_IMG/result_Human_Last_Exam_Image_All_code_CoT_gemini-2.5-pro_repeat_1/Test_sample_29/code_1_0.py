import math

def calculate_hill_number(counts, q):
    """
    Calculates the Hill number for a given list of species counts and an order q.
    """
    total_individuals = sum(counts)
    if total_individuals == 0:
        return 0.0

    proportions = [count / total_individuals for count in counts]

    # Handle the special case for q=1 (Shannon diversity)
    if q == 1:
        shannon_entropy = -sum(p * math.log(p) for p in proportions if p > 0)
        return math.exp(shannon_entropy)
    
    # General formula for q != 1
    else:
        sum_of_proportions_raised_to_q = sum(p**q for p in proportions)
        if sum_of_proportions_raised_to_q == 0:
            return 0.0
        return sum_of_proportions_raised_to_q**(1 / (1 - q))

# Species counts identified from the image
# Species A: 64, Species B: 28, Species C: 24, Species D: 6
species_counts = [64, 28, 24, 6]

# q values as specified in the prompt
q_values = [1, 2, 3, 4]

# Calculate Hill numbers for each q
results = []
for q in q_values:
    hill_val = calculate_hill_number(species_counts, q)
    results.append(hill_val)

# Format the output as four numeric values, rounded to two decimal places, separated by commas.
# The code outputs each number that forms the final answer string.
output_string = f"{results[0]:.2f},{results[1]:.2f},{results[2]:.2f},{results[3]:.2f}"

print(output_string)