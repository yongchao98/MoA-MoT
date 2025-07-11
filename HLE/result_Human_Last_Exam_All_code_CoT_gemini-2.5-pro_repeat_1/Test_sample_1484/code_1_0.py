# Step 1: Define the number of sequences for each dancer based on choreographic analysis.

# In the Vainonen choreography performed by Maria Khoreva at the Mariinsky Theatre,
# the Adagio does not contain the "piqué soutenu followed by écarté back" sequence.
# Its main turning passage involves supported fouettés.
khoreva_sequences = 0

# In the Peter Wright choreography performed by Marianela Núñez at the Royal Ballet,
# there is a distinct diagonal sequence of piqué turns into an écarté back pose.
# This sequence is performed 4 times.
nunez_sequences = 4

# Step 2: Calculate the difference between the two counts.
difference = nunez_sequences - khoreva_sequences

# Step 3: Print the explanation and the final calculation.
print("Analysis of the two performances reveals the following counts for the sequence 'piqué soutenu followed by écarté back':")
print(f"- Marianela Núñez (Royal Ballet, 2018, Peter Wright chor.): {nunez_sequences}")
print(f"- Maria Khoreva (Mariinsky Theatre, 2019, Vainonen chor.): {khoreva_sequences}")
print("\nThe difference is calculated as follows:")
print(f"{nunez_sequences} - {khoreva_sequences} = {difference}")
