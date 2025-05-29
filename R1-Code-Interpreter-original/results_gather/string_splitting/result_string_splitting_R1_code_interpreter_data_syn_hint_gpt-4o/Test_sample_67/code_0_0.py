# Initialize the quantities
machines = {'A': 3, 'B': 2, 'C': 3}
parts = {'X': 0, 'Y': 0, 'Z': 0}

# Define a flag to check if any method was executed in a cycle
method_executed = True

# Continue the cycle until no method can be executed
while method_executed:
    method_executed = False
    
    # Method 1: One unit of machine A -> 2 units of part X and 1 unit of part Y
    if machines['A'] >= 1:
        machines['A'] -= 1
        parts['X'] += 2
        parts['Y'] += 1
        method_executed = True
    
    # Method 2: Two units of machine B -> 1 unit of part X
    if machines['B'] >= 2:
        machines['B'] -= 2
        parts['X'] += 1
        method_executed = True
    
    # Method 3: Two units of machine C -> 1 unit of part Y
    if machines['C'] >= 2:
        machines['C'] -= 2
        parts['Y'] += 1
        method_executed = True
    
    # Method 4: One unit of machine B and one unit of machine C -> 1 unit of machine A
    if machines['B'] >= 1 and machines['C'] >= 1:
        machines['B'] -= 1
        machines['C'] -= 1
        machines['A'] += 1
        method_executed = True
    
    # Method 5: One unit of part X and one unit of part Y -> 1 unit of part Z
    if parts['X'] >= 1 and parts['Y'] >= 1:
        parts['X'] -= 1
        parts['Y'] -= 1
        parts['Z'] += 1
        method_executed = True

# Prepare the final result
result = [str(machines['A']), str(machines['B']), str(machines['C']), str(parts['X']), str(parts['Y']), str(parts['Z'])]
print(result)