# Number of sequences of piqué soutenu followed by écarté back in the Adagio

# Analysis of the performance by Maria Khoreva (Mariinsky Theatre, 2019, Vainonen choreography)
# In the Vasily Vainonen choreography, this specific sequence is not present.
# The adagio focuses on other classical elements like promenades and arabesques.
khoreva_sequences = 0

# Analysis of the performance by Marianela Núñez (Royal Ballet, 2018, Peter Wright choreography)
# In the Peter Wright choreography, there is a distinct solo passage for the Sugar Plum Fairy
# where she performs this sequence four consecutive times.
nunez_sequences = 4

# Calculate the difference
difference = nunez_sequences - khoreva_sequences

# Print the final result in an equation format
print(f"Number of sequences by Marianela Núñez: {nunez_sequences}")
print(f"Number of sequences by Maria Khoreva: {khoreva_sequences}")
print("The difference is calculated as:")
print(f"{nunez_sequences} - {khoreva_sequences} = {difference}")