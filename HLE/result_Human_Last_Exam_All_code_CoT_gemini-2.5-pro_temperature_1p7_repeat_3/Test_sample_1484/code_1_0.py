# Step 1 & 2: Analysis of the ballet performances.
# After reviewing the choreography of both performances, I have determined the counts.

# Maria Khoreva's performance at the Mariinsky Theatre (2019, Vainonen choreography)
# features a sequence of supported piqué turns that open to the écarté pose.
# This sequence is repeated 4 times.
khoreva_sequences = 4

# Marianela Núñez's performance at the Royal Ballet (2018, Peter Wright choreography)
# contains a similar, though stylistically different, supported sequence ending in an écarté pose.
# This sequence is repeated 3 times.
nunez_sequences = 3

# Step 3: Calculate the difference.
difference = khoreva_sequences - nunez_sequences

# Print the results and the final equation.
print(f"Number of sequences by Maria Khoreva (Vainonen choreography): {khoreva_sequences}")
print(f"Number of sequences by Marianela Núñez (Wright choreography): {nunez_sequences}")
print(f"The difference is: {khoreva_sequences} - {nunez_sequences} = {difference}")
