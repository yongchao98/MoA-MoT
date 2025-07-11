# Number of wells with successfully amplified paired heavy and light chains
successful_wells = 4

# Total number of wells into which doublets were sorted
total_sorted_wells = 95

# Calculate the experimental efficiency as a percentage
efficiency = (successful_wells / total_sorted_wells) * 100

# Output the explanation and the calculation
print("The experimental efficiency is calculated by dividing the number of successful outcomes by the total number of attempts.")
print(f"In this case, the equation is: ({successful_wells} / {total_sorted_wells}) * 100")
print(f"The resulting efficiency is approximately {efficiency:.2f}%.")