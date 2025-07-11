# The Nutcracker Adagio Choreography Analysis

# Number of repeating piqué balance sequences in the Vasily Vainonen
# choreography, as performed by Maria Khoreva at the Mariinsky Theatre.
khoreva_sequences = 4

# Number of repeating piqué balance sequences in the Peter Wright
# choreography, as performed by Marianela Núñez at the Royal Ballet.
nunez_sequences = 4

# Calculate the difference between the two counts.
difference = khoreva_sequences - nunez_sequences

# Display the final calculation.
print(f"Maria Khoreva's count: {khoreva_sequences}")
print(f"Marianela Núñez's count: {nunez_sequences}")
print("The difference is calculated as:")
print(f"{khoreva_sequences} - {nunez_sequences} = {difference}")