# This code summarizes the logical deduction of the substituents based on the
# known backbone rearrangement mechanism for triterpenoids.

# Dictionary to store the identified substituents.
# V, W, X, Y, Z corresponding to positions 1, 2, 3, 4, 5.
substituents = {
    '1': 'CH3', # Result for position 1 (at C-4)
    '2': 'H',   # Result for position 2 (at C-10)
    '3': 'CH3', # Result for position 3 (at C-8)
    '4': 'H',   # Result for position 4 (at C-14)
    '5': 'H'    # Result for position 5 (at C-18)
}

# Print the final answer in the required format.
print(f"1 = {substituents['1']}, 2 = {substituents['2']}, 3 = {substituents['3']}, 4 = {substituents['4']}, 5 = {substituents['5']}")