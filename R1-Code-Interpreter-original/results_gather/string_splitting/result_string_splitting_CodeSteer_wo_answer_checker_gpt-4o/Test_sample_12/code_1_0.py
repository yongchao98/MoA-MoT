# Initial quantities
machines = {'A': 1, 'B': 1, 'C': 3}
parts = {'X': 0, 'Y': 0, 'Z': 0}

# Function to execute the methods
def execute_methods():
    while True:
        executed = False
        
        # Method 1
        if machines['A'] >= 1:
            machines['A'] -= 1
            parts['X'] += 2
            parts['Y'] += 1
            executed = True
        
        # Method 2
        elif machines['B'] >= 2:
            machines['B'] -= 2
            parts['X'] += 1
            executed = True
        
        # Method 3
        elif machines['C'] >= 2:
            machines['C'] -= 2
            parts['Y'] += 1
            executed = True
        
        # Method 4
        elif machines['B'] >= 1 and machines['C'] >= 1:
            machines['B'] -= 1
            machines['C'] -= 1
            machines['A'] += 1
            executed = True
        
        # Method 5
        elif parts['X'] >= 1 and parts['Y'] >= 1:
            parts['X'] -= 1
            parts['Y'] -= 1
            parts['Z'] += 1
            executed = True
        
        # If no method was executed, break the loop
        if not executed:
            break

# Execute the methods
execute_methods()

# Prepare the result
result = [str(machines['A']), str(machines['B']), str(machines['C']),
          str(parts['X']), str(parts['Y']), str(parts['Z'])]

# Print the result
print(result)