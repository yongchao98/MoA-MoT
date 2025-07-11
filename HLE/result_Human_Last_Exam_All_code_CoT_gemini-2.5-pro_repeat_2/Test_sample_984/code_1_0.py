# Define the number of successful outcomes and total trials
successful_wells = 4
total_sorted_wells = 95

# Calculate the efficiency
efficiency = (successful_wells / total_sorted_wells) * 100

# Print the final result in a clear sentence, showing the numbers used
print(f"The experiment yielded paired heavy and light chains in {successful_wells} out of {total_sorted_wells} sorted wells.")
print(f"The equation for the efficiency is: ({successful_wells} / {total_sorted_wells}) * 100")
print(f"The experimental efficiency is: {efficiency:.2f}%")