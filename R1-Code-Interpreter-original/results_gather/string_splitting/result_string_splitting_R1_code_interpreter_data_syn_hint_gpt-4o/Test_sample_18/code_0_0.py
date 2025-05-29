# Initial counts
machines = {'A': 1, 'B': 2, 'C': 3}
parts = {'X': 0, 'Y': 0, 'Z': 0}

# Function to simulate the process
def simulate_process(machines, parts):
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

simulate_process(machines, parts)

# Output the final counts
result = [str(machines['A']), str(machines['B']), str(machines['C']), str(parts['X']), str(parts['Y']), str(parts['Z'])]
print(result)