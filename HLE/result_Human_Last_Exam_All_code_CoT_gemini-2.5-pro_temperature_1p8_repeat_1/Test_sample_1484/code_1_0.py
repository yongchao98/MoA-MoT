# Number of "piqué soutenu followed by écarté back" sequences for each dancer.

# Maria Khoreva in the Vainonen choreography (2019, Mariinsky Theatre)
# This choreography does not include the specified sequence.
khoreva_sequences = 0

# Marianela Núñez in the Peter Wright choreography (2018, Royal Ballet)
# This choreography includes a diagonal where this sequence is performed four times.
nunez_sequences = 4

# Calculate the difference
difference = nunez_sequences - khoreva_sequences

# Print the final result showing the numbers used in the equation
print(f"To find the difference, we subtract the number of sequences performed by Maria Khoreva from that of Marianela Núñez.")
print(f"Calculation: {nunez_sequences} - {khoreva_sequences} = {difference}")