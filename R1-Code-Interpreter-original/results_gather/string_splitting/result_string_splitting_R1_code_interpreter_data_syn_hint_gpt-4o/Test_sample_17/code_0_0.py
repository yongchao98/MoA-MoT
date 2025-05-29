# Initialize the counts of machines and parts
machines = {'A': 3, 'B': 3, 'C': 3}
parts = {'X': 0, 'Y': 0, 'Z': 0}

# Function to apply the methods in sequence
def apply_methods(machines, parts):
    while True:
        # Method 1
        if machines['A'] >= 1:
            machines['A'] -= 1
            parts['X'] += 2
            parts['Y'] += 1
        # Method 2
        elif machines['B'] >= 2:
            machines['B'] -= 2
            parts['X'] += 1
        # Method 3
        elif machines['C'] >= 2:
            machines['C'] -= 2
            parts['Y'] += 1
        # Method 4
        elif machines['B'] >= 1 and machines['C'] >= 1:
            machines['B'] -= 1
            machines['C'] -= 1
            machines['A'] += 1
        # Method 5
        elif parts['X'] >= 1 and parts['Y'] >= 1:
            parts['X'] -= 1
            parts['Y'] -= 1
            parts['Z'] += 1
        else:
            break

# Apply the methods
apply_methods(machines, parts)

# Prepare the result in the required format
result = [str(machines['A']), str(machines['B']), str(machines['C']),
          str(parts['X']), str(parts['Y']), str(parts['Z'])]

# Print the result
print(result)