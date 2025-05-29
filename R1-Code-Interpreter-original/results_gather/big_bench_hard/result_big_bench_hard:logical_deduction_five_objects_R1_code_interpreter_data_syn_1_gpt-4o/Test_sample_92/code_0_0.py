# Define the positions
positions = ['Ana', 'Amy', 'Eli', 'Mya', 'Dan']

# Check the conditions
condition1 = positions.index('Amy') < positions.index('Eli')  # Amy finished above Eli
condition2 = positions.index('Ana') < positions.index('Amy')  # Ana finished above Amy
condition3 = positions.index('Mya') > positions.index('Eli')  # Mya finished below Eli
condition4 = positions.index('Dan') == 4                      # Dan finished last

# Verify all conditions
all_conditions_met = condition1 and condition2 and condition3 and condition4

# Output the result
print("All conditions met:", all_conditions_met)
print("Eli's position:", positions.index('Eli') + 1)