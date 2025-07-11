# Number of "piqué soutenu followed by écarté back" sequences in each performance

# In the 2019 Mariinsky performance choreographed by Vasily Vainonen,
# Maria Khoreva does not perform this specific sequence.
khoreva_sequences = 0

# In the 2018 Royal Ballet performance choreographed by Peter Wright,
# Marianela Núñez famously performs this sequence 4 times in a diagonal.
nunez_sequences = 4

# Calculate the difference
difference = nunez_sequences - khoreva_sequences

# Print the result in a full sentence showing the equation
print(f"Number of sequences by Marianela Núñez: {nunez_sequences}")
print(f"Number of sequences by Maria Khoreva: {khoreva_sequences}")
print(f"The difference is {nunez_sequences} - {khoreva_sequences} = {difference}")
