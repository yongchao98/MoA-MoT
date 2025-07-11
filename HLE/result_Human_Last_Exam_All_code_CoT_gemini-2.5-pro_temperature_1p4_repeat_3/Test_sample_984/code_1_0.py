import sys

# The number of wells where paired heavy and light chains were successfully amplified.
successful_wells = 4

# The total number of wells into which doublets were sorted (excluding the negative control).
total_sorted_wells = 95

# The efficiency is the ratio of successful outcomes to the total number of attempts, expressed as a percentage.
# We check if total_sorted_wells is zero to avoid division by zero error, although it's not the case here.
if total_sorted_wells > 0:
    efficiency_percentage = (successful_wells / total_sorted_wells) * 100
else:
    efficiency_percentage = 0
    # Exit if the denominator is 0, as the calculation would be invalid.
    print("Error: Total number of sorted wells cannot be zero.", file=sys.stderr)
    sys.exit(1)

# To fulfill the request, we will print out each number and the final equation.
print(f"Number of successful wells: {successful_wells}")
print(f"Total number of sorted wells: {total_sorted_wells}")
print(f"Efficiency calculation: ({successful_wells} / {total_sorted_wells}) * 100%")
print(f"The final efficiency is: {efficiency_percentage:.2f}%")