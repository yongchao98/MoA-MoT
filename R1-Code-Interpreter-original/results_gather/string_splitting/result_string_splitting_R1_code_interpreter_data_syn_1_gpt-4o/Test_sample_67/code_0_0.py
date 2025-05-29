# Initial quantities
machines = {'A': 3, 'B': 2, 'C': 3}
parts = {'X': 0, 'Y': 0, 'Z': 0}

# Function to apply methods
def apply_methods(machines, parts):
    while True:
        # Method 1: One unit of machine A -> 2 units of part X and 1 unit of part Y
        if machines['A'] >= 1:
            machines['A'] -= 1
            parts['X'] += 2
            parts['Y'] += 1
        # Method 2: Two units of machine B -> 1 unit of part X
        elif machines['B'] >= 2:
            machines['B'] -= 2
            parts['X'] += 1
        # Method 3: Two units of machine C -> 1 unit of part Y
        elif machines['C'] >= 2:
            machines['C'] -= 2
            parts['Y'] += 1
        # Method 4: One unit of machine B and one unit of machine C -> 1 unit of machine A
        elif machines['B'] >= 1 and machines['C'] >= 1:
            machines['B'] -= 1
            machines['C'] -= 1
            machines['A'] += 1
        # Method 5: One unit of part X and one unit of part Y -> 1 unit of part Z
        elif parts['X'] >= 1 and parts['Y'] >= 1:
            parts['X'] -= 1
            parts['Y'] -= 1
            parts['Z'] += 1
        else:
            # No more methods can be executed
            break

apply_methods(machines, parts)

# Output the final quantities
result = [str(machines['A']), str(machines['B']), str(machines['C']), str(parts['X']), str(parts['Y']), str(parts['Z'])]
print(result)