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
    machines = {'A': 4, 'B': 3, 'C': 4}
    parts = {'X': 0, 'Y': 0, 'Z': 0}
    current_method = 1
    
    while True:
        executed_any = False
        for i in range(5):
            current_method = i + 1
            while can_execute_method(current_method, machines, parts):
                execute_method(current_method, machines, parts)
                executed_any = True
                print(f"After Method {current_method}:", 
                      f"Machines: A={machines['A']}, B={machines['B']}, C={machines['C']}", 
                      f"Parts: X={parts['X']}, Y={parts['Y']}, Z={parts['Z']}")
        
        if not executed_any:
            break
    
    result = [str(machines['A']), str(machines['B']), str(machines['C']), 
             str(parts['X']), str(parts['Y']), str(parts['Z'])]
    print("\nFinal result:", result)

simulate_dismantling()