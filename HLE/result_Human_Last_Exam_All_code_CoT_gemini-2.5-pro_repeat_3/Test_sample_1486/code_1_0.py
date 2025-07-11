# Plan:
# 1. Identify the number of cabrioles devants in the first diagonal sequence.
# 2. Identify the number of cabrioles devants in the second diagonal sequence.
# 3. Sum these numbers to get the total.
# 4. Print the final equation and the result.

# In Maria Khoreva's solo, she performs two identical diagonal passes of jumps.
# Each pass includes two cabrioles devants.
cabrioles_in_first_diagonal = 2
cabrioles_in_second_diagonal = 2

# Calculate the total number of cabrioles
total_cabrioles = cabrioles_in_first_diagonal + cabrioles_in_second_diagonal

# Print the equation showing the breakdown and the final count.
print(f"Maria Khoreva executed {cabrioles_in_first_diagonal} cabrioles in the first diagonal sequence and {cabrioles_in_second_diagonal} in the second.")
print("The total number of cabrioles devants is calculated as follows:")
print(f"{cabrioles_in_first_diagonal} + {cabrioles_in_second_diagonal} = {total_cabrioles}")