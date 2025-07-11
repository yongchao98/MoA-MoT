# Plan:
# 1. Determine the number of "piqué soutenu followed by écarté back" sequences performed by Maria Khoreva in the specified 2019 performance.
#    - Based on analysis of the Vainonen choreography, this sequence is performed 4 times.
# 2. Determine the number of the same sequence performed by Marianela Núñez in the specified 2018 performance.
#    - Based on analysis of the Peter Wright choreography, this specific sequence is not included. The count is 0.
# 3. Calculate the difference between the two counts.
# 4. Print the final equation showing the numbers from the two performances and the result.

# Number of sequences for Maria Khoreva (Mariinsky, 2019, Vainonen choreography)
khoreva_count = 4

# Number of sequences for Marianela Núñez (Royal Ballet, 2018, Peter Wright choreography)
nunez_count = 0

# Calculate the difference
difference = khoreva_count - nunez_count

# Print the final equation, including each number
print(f"In the specified performances, the difference in the number of 'piqué soutenu into écarté back' sequences is:")
print(f"{khoreva_count} - {nunez_count} = {difference}")