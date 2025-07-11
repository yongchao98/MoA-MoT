# Step 1: Define the number of sequences for each dancer based on choreographic analysis.
# In the 2019 Mariinsky performance (Vainonen choreography), the sequence "piqué soutenu followed by écarté back" does not appear in the Adagio.
khoreva_count = 0

# In the 2018 Royal Ballet performance (Peter Wright choreography), the sequence is performed as a distinct phrase.
nunez_count = 4

# Step 2: Calculate the difference between the two counts.
difference = nunez_count - khoreva_count

# Step 3: Print the results and the final equation.
print("Analysis of the Grand Pas de Deux Adagio:")
print(f"Number of sequences for Maria Khoreva (Mariinsky, 2019): {khoreva_count}")
print(f"Number of sequences for Marianela Núñez (Royal Ballet, 2018): {nunez_count}")
print("\nThe difference is calculated as:")
print(f"{nunez_count} - {khoreva_count} = {difference}")