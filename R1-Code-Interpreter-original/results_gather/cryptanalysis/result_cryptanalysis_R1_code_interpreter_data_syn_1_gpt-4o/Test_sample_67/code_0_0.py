# Define the possible numbers and letters based on the conditions
numbers = [4, 5]
letters = ['L', 'M']

# Check each condition
conditions_met = [
    # Guess 1: 79FV
    (79, 'FV', False),
    # Guess 2: 32PZ
    (32, 'PZ', False),
    # Guess 3: 09EF
    (9, 'EF', False),
    # Guess 4: 58QD
    (58, 'QD', False),
    # Guess 5: 79OJ
    (79, 'OJ', False),
    # Guess 6: 64LY
    (64, 'LY', False),
    # Guess 7: 48HI
    (48, 'HI', False),
    # Guess 8: 45TG
    (45, 'TG', True),
    # Guess 9: 31IB
    (31, 'IB', False),
    # Guess 10: 94VW
    (94, 'VW', False),
    # Guess 11: 70XN
    (70, 'XN', False),
    # Guess 12: 70BI
    (70, 'BI', False),
    # Guess 13: 89UG
    (89, 'UG', False),
    # Guess 14: 70KG
    (70, 'KG', False)
]

# Check if the combination 45LM meets all conditions
combination = (45, 'LM')
is_correct = all(
    (combination[0] != guess[0] or combination[1] != guess[1]) == guess[2]
    for guess in conditions_met
)

print(is_correct)