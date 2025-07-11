from fractions import Fraction

def calculate_moment(k):
    """
    Calculates the k-th moment of the dimensionless conductance g
    for a disordered Majorana wire at the critical point using the formula <g^k> = 1 / (2k + 1).
    """
    return Fraction(1, 2 * k + 1)

# Define the order of the moments we need to calculate
k_avg = 1       # For the average value <g>
k_fourth = 4    # For the fourth moment <g^4>

# Calculate the average conductance <g>
avg_g = calculate_moment(k_avg)

# Calculate the fourth statistical moment <g^4>
fourth_moment_g = calculate_moment(k_fourth)

# Calculate the ratio
ratio = fourth_moment_g / avg_g

print("The problem asks for the ratio of the fourth moment of conductance to its average value.")
print("For a disordered Majorana wire at the critical point, the k-th moment <g^k> follows the universal formula: 1 / (2k + 1).")
print("-" * 70)

# Print the calculation for the average value
print(f"Step 1: Calculate the average conductance, <g> (the first moment, k={k_avg}).")
print(f"   <g> = 1 / (2 * {k_avg} + 1) = {avg_g}")
print()

# Print the calculation for the fourth moment
print(f"Step 2: Calculate the fourth moment, <g^4> (k={k_fourth}).")
print(f"   <g^4> = 1 / (2 * {k_fourth} + 1) = {fourth_moment_g}")
print()

# Print the final ratio calculation
print("Step 3: Calculate the final ratio <g^4> / <g>.")
print(f"   Ratio = <g^4> / <g> = ({fourth_moment_g}) / ({avg_g}) = {ratio}")
print()
print(f"The final result for the ratio is {ratio}.")
