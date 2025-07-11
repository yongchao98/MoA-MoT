# Step 1: Define the number of sequences for each performance based on choreographic analysis.

# In the Adagio of the 2018 Royal Ballet performance (Peter Wright choreography),
# Marianela Núñez performs the sequence of "piqué soutenu followed by écarté back" four times.
nunez_sequences = 4

# In the Adagio of the 2019 Mariinsky Theatre performance (Vasily Vainonen choreography),
# Maria Khoreva does not perform this specific sequence.
khoreva_sequences = 0

# Step 2: Calculate the difference between the two counts.
difference = nunez_sequences - khoreva_sequences

# Step 3: Print the result clearly, showing each number in the final equation.
print(f"Number of sequences by Marianela Núñez: {nunez_sequences}")
print(f"Number of sequences by Maria Khoreva: {khoreva_sequences}")
print(f"The difference is: {nunez_sequences} - {khoreva_sequences} = {difference}")