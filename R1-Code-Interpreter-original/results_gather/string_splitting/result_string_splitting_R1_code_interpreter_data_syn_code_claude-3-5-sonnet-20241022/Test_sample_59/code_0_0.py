def can_apply_method(method, machines, parts):
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

def apply_method(method, machines, parts):
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
    machines = {'A': 4, 'B': 4, 'C': 4}
    parts = {'X': 0, 'Y': 0, 'Z': 0}
    
    while True:
        any_method_applied = False
        for method in range(1, 6):
            if can_apply_method(method, machines, parts):
                apply_method(method, machines, parts)
                any_method_applied = True
        
        if not any_method_applied:
            break
    
    result = [str(machines['A']), str(machines['B']), str(machines['C']), 
              str(parts['X']), str(parts['Y']), str(parts['Z'])]
    print(result)

simulate_dismantling()