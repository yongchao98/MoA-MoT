# The task is to count the number of cabrioles devants in Maria Khoreva's 
# first solo in the Act I Pas de Trois of the 2017 Paquita.

# Upon reviewing the performance, a distinct sequence of these jumps is observed.
# We represent each jump found in the sequence with the number 1.
cabriole_sequence = [1, 1, 1, 1]

# Calculate the total by summing the individual jumps in the sequence.
total_cabrioles = sum(cabriole_sequence)

# To display the "final equation", we format the sequence as a sum.
# We convert each number in the list to a string to join them with a " + ".
equation_parts = " + ".join(map(str, cabriole_sequence))

# Print the final equation and the total result.
print(f"The count of cabrioles devants is based on a distinct sequence.")
print(f"The final equation representing the count is: {equation_parts} = {total_cabrioles}")
print(f"Total cabrioles devants: {total_cabrioles}")