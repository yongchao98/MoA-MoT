# Initialize the quantities of machines and parts
machines = {'A': 3, 'B': 5, 'C': 6}
parts = {'X': 0, 'Y': 0, 'Z': 0}

# Define a flag to check if any method was executed in the last cycle
method_executed = True

# Continue cycling through the methods until no method can be executed
while method_executed:
    method_executed = False
    
    # Method 1: One unit of machine A can be dismantled into two units of part X and one unit of part Y
    if machines['A'] >= 1:
        machines['A'] -= 1
        parts['X'] += 2
        parts['Y'] += 1
        method_executed = True
    
    # Method 2: Two units of machine B can be dismantled into one unit of part X
    elif machines['B'] >= 2:
        machines['B'] -= 2
        parts['X'] += 1
        method_executed = True
    
    # Method 3: Dismantling two machine C units yields one part Y
    elif machines['C'] >= 2:
        machines['C'] -= 2
        parts['Y'] += 1
        method_executed = True
    
    # Method 4: One unit of machine B and one unit of machine C can be combined into one unit of machine A
    elif machines['B'] >= 1 and machines['C'] >= 1:
        machines['B'] -= 1
        machines['C'] -= 1
        machines['A'] += 1
        method_executed = True
    
    # Method 5: One unit of part X and one unit of part Y can be combined into one unit of part Z
    elif parts['X'] >= 1 and parts['Y'] >= 1:
        parts['X'] -= 1
        parts['Y'] -= 1
        parts['Z'] += 1
        method_executed = True

# Prepare the final result in the specified format
result = [str(machines['A']), str(machines['B']), str(machines['C']),
          str(parts['X']), str(parts['Y']), str(parts['Z'])]

# Print the result
print(result)