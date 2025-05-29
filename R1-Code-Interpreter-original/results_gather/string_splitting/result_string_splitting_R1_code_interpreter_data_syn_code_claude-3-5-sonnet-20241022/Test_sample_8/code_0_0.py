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
    if method == 1:
        machines['A'] -= 1
        parts['X'] += 2
        parts['Y'] += 1
    elif method == 2:
        machines['B'] -= 2
        parts['X'] += 1
    elif method == 3:
        machines['C'] -= 2
        parts['Y'] += 1
    elif method == 4:
        machines['B'] -= 1
        machines['C'] -= 1
        machines['A'] += 1
    elif method == 5:
        parts['X'] -= 1
        parts['Y'] -= 1
        parts['Z'] += 1

def simulate_dismantling():
    machines = {'A': 2, 'B': 2, 'C': 3}
    parts = {'X': 0, 'Y': 0, 'Z': 0}
    current_method = 1
    made_progress = True
    
    while made_progress:
        made_progress = False
        for _ in range(5):  # Try each method once per cycle
            if can_execute_method(current_method, machines, parts):
                execute_method(current_method, machines, parts)
                made_progress = True
            current_method = (current_method % 5) + 1
    
    result = [str(machines['A']), str(machines['B']), str(machines['C']), 
             str(parts['X']), str(parts['Y']), str(parts['Z'])]
    print(result)

simulate_dismantling()