# Number of wells sorted based on the doublet gate
total_wells_sorted = 95

# Number of wells that yielded paired heavy and light chains after PCR
successful_wells = 4

# Calculate the efficiency
# The formula is (successful events / total events) * 100
efficiency_percentage = (successful_wells / total_wells_sorted) * 100

# Print the result in a clear format, showing the numbers used in the equation
print(f"The observed experimental efficiency is calculated from the number of successful wells and the total number of sorted wells.")
print(f"Equation: {successful_wells} / {total_wells_sorted} * 100")
print(f"Result: {efficiency_percentage:.2f}%")