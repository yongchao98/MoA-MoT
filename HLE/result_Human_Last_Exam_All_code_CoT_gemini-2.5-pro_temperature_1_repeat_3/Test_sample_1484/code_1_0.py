# Number of sequences performed by Maria Khoreva in the Vainonen choreography.
# This choreography does not contain the specified "piqué soutenu followed by écarté back" sequence.
khoreva_sequences = 0

# Number of sequences performed by Marianela Núñez in the Peter Wright choreography.
# This choreography contains a similar step, but the pose is écarté devant (front), not écarté back.
nunez_sequences = 0

# Calculate the difference between the two counts.
difference = abs(khoreva_sequences - nunez_sequences)

# Print the final equation showing each number.
print(f"Maria Khoreva's count: {khoreva_sequences}")
print(f"Marianela Núñez's count: {nunez_sequences}")
print(f"The difference is: {khoreva_sequences} - {nunez_sequences} = {difference}")