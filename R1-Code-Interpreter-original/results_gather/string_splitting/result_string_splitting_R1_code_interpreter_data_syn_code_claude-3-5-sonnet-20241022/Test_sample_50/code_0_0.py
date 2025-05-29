def can_execute_method(method, machines, parts):
    if method == 1:
        return machines['A'] >= 1
    elif method == 2:
        return machines['B'] >= 2
    elif method == 3:
        return machines['C'] >= 2
    elif method == 4:
        return machines['B'] >= 1 and machines['C'] >= 1
    elif method == 5:
        return parts['X'] >= 1 and parts['Y'] >= 1
    return False

def execute_method(method, machines, parts):
    if method == 1:  # A -> 2X + Y
        machines['A'] -= 1
        parts['X'] += 2
        parts['Y'] += 1
    elif method == 2:  # 2B -> X
        machines['B'] -= 2
        parts['X'] += 1
    elif method == 3:  # 2C -> Y
        machines['C'] -= 2
        parts['Y'] += 1
    elif method == 4:  # B + C -> A
        machines['B'] -= 1
        machines['C'] -= 1
        machines['A'] += 1
    elif method == 5:  # X + Y -> Z
        parts['X'] -= 1
        parts['Y'] -= 1
        parts['Z'] += 1

# Initialize quantities
machines = {'A': 2, 'B': 4, 'C': 3}
parts = {'X': 0, 'Y': 0, 'Z': 0}

# Continue until no method can be executed
while True:
    any_method_executed = False
    for method in range(1, 6):
        if can_execute_method(method, machines, parts):
            execute_method(method, machines, parts)
            any_method_executed = True
    if not any_method_executed:
        break

# Format the result
result = [str(machines['A']), str(machines['B']), str(machines['C']), 
          str(parts['X']), str(parts['Y']), str(parts['Z'])]
print(result)