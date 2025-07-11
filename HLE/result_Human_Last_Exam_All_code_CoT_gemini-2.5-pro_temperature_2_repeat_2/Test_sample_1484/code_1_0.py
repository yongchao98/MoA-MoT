# The number of times the sequence is performed in each ballet version.
# In Vasily Vainonen's choreography performed by Maria Khoreva, this specific repeating sequence is not present.
khoreva_vainonen_count = 0

# In Peter Wright's choreography performed by Marianela Núñez, this sequence is famously repeated on a diagonal.
nunez_wright_count = 4

# Calculate the difference.
difference = nunez_wright_count - khoreva_vainonen_count

# Print the final result, showing the equation.
print(f"The number of 'piqué soutenu into écarté back' sequences performed by Marianela Núñez (Wright choreography) is {nunez_wright_count}.")
print(f"The number of 'piqué soutenu into écarté back' sequences performed by Maria Khoreva (Vainonen choreography) is {khoreva_vainonen_count}.")
print(f"The difference is {nunez_wright_count} - {khoreva_vainonen_count} = {difference}.")