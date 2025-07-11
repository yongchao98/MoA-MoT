# Define the relationship between outer radius (ro) and inner radius (ri)
# as given in the problem statement.
# ro = 2 * ri
# To simplify the calculation, we can use normalized values where ri=1, so ro=2.
ri = 1
ro = 2

# The formula for the ratio of the maximum tangential stress (σ_max) to the 
# internal pressure (P_i) in a thick-walled cylinder is:
# Ratio = (ro^2 + ri^2) / (ro^2 - ri^2)

# Calculate the squared terms
ro_squared = ro**2
ri_squared = ri**2

# Calculate the numerator and the denominator of the ratio formula
numerator = ro_squared + ri_squared
denominator = ro_squared - ri_squared

# Display the calculation process, showing each number in the final equation
print("The ratio of maximum tangential stress (σ_max) to internal pressure (P_i) is:")
print("σ_max / P_i = (ro² + ri²) / (ro² - ri²)")
print("\nGiven that ro = 2 * ri, we can set ro=2 and ri=1 for the calculation.")
print(f"σ_max / P_i = ({ro}² + {ri}²) / ({ro}² - {ri}²)")
print(f"σ_max / P_i = ({ro_squared} + {ri_squared}) / ({ro_squared} - {ri_squared})")
print(f"σ_max / P_i = {int(numerator)} / {int(denominator)}")

# The final ratio as a decimal
final_ratio_decimal = numerator / denominator
print(f"\nThe resulting ratio is 5/3, which is approximately {final_ratio_decimal:.3f}.")
