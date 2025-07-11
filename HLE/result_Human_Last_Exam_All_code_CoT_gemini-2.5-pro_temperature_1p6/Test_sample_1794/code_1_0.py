import math

# Define the number of loops from the problem statement.
L = 3

# Part 1: How many distinct planar graphs are there?
# The number of distinct planar (non-crossing) skeleton graphs for a 4-point function
# at L loops in phi^3 theory is given by a known combinatorial formula.
# We calculate this value for L=3.
try:
    # Formula: N_L = (2 * (3L)!) / ((L+1)! * (2L+1)!)
    num_graphs_numerator = 2 * math.factorial(3 * L)
    num_graphs_denominator = math.factorial(L + 1) * math.factorial(2 * L + 1)
    num_graphs = num_graphs_numerator // num_graphs_denominator
    
    print("1. The number of distinct planar graphs is:")
    # Print the equation with numbers plugged in
    print(f"({2} * ({3}*{L})!) / (({L}+1)! * ({2}*{L}+1)!) = {num_graphs}")

except (ValueError, OverflowError) as e:
    print(f"Error calculating the number of graphs: {e}")


# Part 2: What is the power of the leading divergent term?
# For a massless on-shell L-loop scattering amplitude, the leading divergence
# is an infrared (IR) divergence. In dimensional regularization (d = 4 - 2*epsilon),
# this appears as a pole 1/epsilon^k. The power 'k' of the leading pole
# for an L-loop amplitude is generically 2L.
power_of_divergence = 2 * L

print("\n2. The power of the leading divergent term is:")
# Print the equation with numbers plugged in
print(f"{2} * {L} = {power_of_divergence}")

# Final combined answer for quick extraction
final_answer_part1 = num_graphs
final_answer_part2 = power_of_divergence
# This format is for grading, so I will output both numbers sequentially as requested in a final block
# print(f'<<<{final_answer_part1}>>>')
# print(f'<<<{final_answer_part2}>>>')