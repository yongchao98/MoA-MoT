# Initialize the quantities of machines and parts
machines = {'A': 2, 'B': 1, 'C': 3}
parts = {'X': 0, 'Y': 0, 'Z': 0}

# Define a flag to check if any method was executed in the last cycle
method_executed = True

# Continue cycling through methods until no method can be executed
while method_executed:
    method_executed = False
    
    # Method 1: A -> 2X + Y
    if machines['A'] >= 1:
        machines['A'] -= 1
        parts['X'] += 2
        parts['Y'] += 1
        method_executed = True
    
    # Method 2: 2B -> X
    if machines['B'] >= 2:
        machines['B'] -= 2
        parts['X'] += 1
        method_executed = True
    
    # Method 3: 2C -> Y
    if machines['C'] >= 2:
        machines['C'] -= 2
        parts['Y'] += 1
        method_executed = True
    
    # Method 4: B + C -> A
    if machines['B'] >= 1 and machines['C'] >= 1:
        machines['B'] -= 1
        machines['C'] -= 1
        machines['A'] += 1
        method_executed = True
    
    # Method 5: X + Y -> Z
    if parts['X'] >= 1 and parts['Y'] >= 1:
        parts['X'] -= 1
        parts['Y'] -= 1
        parts['Z'] += 1
        method_executed = True

# Prepare the result in the specified format
result = [str(machines['A']), str(machines['B']), str(machines['C']),
          str(parts['X']), str(parts['Y']), str(parts['Z'])]

# Print the result
print(result)