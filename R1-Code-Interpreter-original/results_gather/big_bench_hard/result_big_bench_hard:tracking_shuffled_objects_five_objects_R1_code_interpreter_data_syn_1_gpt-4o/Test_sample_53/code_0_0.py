# Initial positions
positions = {
    'Alice': 'Center Midfielder',
    'Bob': 'Fullback',
    'Claire': 'Goalkeeper',
    'Dave': 'Left Winger',
    'Eve': 'Right Winger'
}

# Swap 1: Eve and Claire trade positions
positions['Eve'], positions['Claire'] = positions['Claire'], positions['Eve']

# Swap 2: Alice and Dave trade positions
positions['Alice'], positions['Dave'] = positions['Dave'], positions['Alice']

# Swap 3: Claire and Bob trade positions
positions['Claire'], positions['Bob'] = positions['Bob'], positions['Claire']

# Swap 4: Eve and Alice trade positions
positions['Eve'], positions['Alice'] = positions['Alice'], positions['Eve']

# Swap 5: Claire and Eve trade positions
positions['Claire'], positions['Eve'] = positions['Eve'], positions['Claire']

# Final position of Eve
print(positions['Eve'])