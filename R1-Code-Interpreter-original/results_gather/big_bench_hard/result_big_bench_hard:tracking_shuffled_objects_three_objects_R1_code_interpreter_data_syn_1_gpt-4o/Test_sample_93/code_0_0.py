# Initial positions
positions = {
    'Alice': 'Right Winger',
    'Bob': 'Left Midfielder',
    'Claire': 'Center Midfielder'
}

# First swap: Claire and Alice
positions['Alice'], positions['Claire'] = positions['Claire'], positions['Alice']

# Second swap: Bob and Alice
positions['Alice'], positions['Bob'] = positions['Bob'], positions['Alice']

# Third swap: Claire and Bob
positions['Claire'], positions['Bob'] = positions['Bob'], positions['Claire']

# Final position of Alice
print(positions['Alice'])