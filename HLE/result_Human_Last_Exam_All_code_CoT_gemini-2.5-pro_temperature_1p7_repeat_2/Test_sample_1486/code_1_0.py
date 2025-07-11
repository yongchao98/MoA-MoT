# In Maria Khoreva's Act I pas de trois variation from the 2017 "Paquita",
# she performs a specific diagonal sequence of jumps multiple times.

# Number of times the sequence is performed
num_sequences = 4

# Number of cabrioles devants executed in each sequence
cabrioles_per_sequence = 2

# To find the total, we multiply the number of sequences
# by the number of cabrioles performed in each one.
total_cabrioles = num_sequences * cabrioles_per_sequence

# The final equation is:
print(f"{num_sequences} sequences * {cabrioles_per_sequence} cabrioles per sequence = {total_cabrioles} total cabrioles devants")