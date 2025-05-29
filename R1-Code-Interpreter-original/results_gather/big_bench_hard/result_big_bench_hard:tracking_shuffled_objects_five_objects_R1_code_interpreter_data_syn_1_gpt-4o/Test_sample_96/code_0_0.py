# Initial ball positions
positions = {
    'Alice': 'Brown',
    'Bob': 'Blue',
    'Claire': 'Orange',
    'Dave': 'Pink',
    'Eve': 'White'
}

# Swap operations
swaps = [
    ('Dave', 'Alice'),
    ('Alice', 'Claire'),
    ('Claire', 'Bob'),
    ('Eve', 'Bob'),
    ('Alice', 'Dave')
]

# Perform swaps
for player1, player2 in swaps:
    positions[player1], positions[player2] = positions[player2], positions[player1]

# Output the ball Claire has at the end
print(positions['Claire'])