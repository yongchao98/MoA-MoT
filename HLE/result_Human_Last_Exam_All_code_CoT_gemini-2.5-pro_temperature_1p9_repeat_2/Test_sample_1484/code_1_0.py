# Number of "piqué soutenu followed by écarté back" sequences

# In Marianela Núñez's 2018 performance (Peter Wright choreography),
# this sequence is a signature part of the Adagio, performed four times.
nunez_sequences = 4

# In Maria Khoreva's 2019 performance (Vasily Vainonen choreography),
# this specific sequence is not included in the Adagio.
khoreva_sequences = 0

# Calculate the difference
difference = nunez_sequences - khoreva_sequences

print("Analyzing the Adagio from the Sugar Plum grand pas de deux:")
print(f"Number of sequences by Marianela Núñez (2018, Royal Ballet): {nunez_sequences}")
print(f"Number of sequences by Maria Khoreva (2019, Mariinsky Theatre): {khoreva_sequences}")
print("\nCalculating the difference:")
print(f"{nunez_sequences} - {khoreva_sequences} = {difference}")
