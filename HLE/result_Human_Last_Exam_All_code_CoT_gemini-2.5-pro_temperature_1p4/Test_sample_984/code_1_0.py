# Define the number of successful outcomes and total trials from the experiment
successful_wells = 4
total_sorted_wells = 95

# Calculate the efficiency of the experiment as a percentage
efficiency_percentage = (successful_wells / total_sorted_wells) * 100

# Print the calculation and the result in a clear, descriptive sentence.
# The f-string formats the floating-point number to two decimal places.
print(f"The experimental efficiency was extremely low. The calculation is:")
print(f"({successful_wells} successful wells / {total_sorted_wells} sorted wells) * 100 = {efficiency_percentage:.2f}%")