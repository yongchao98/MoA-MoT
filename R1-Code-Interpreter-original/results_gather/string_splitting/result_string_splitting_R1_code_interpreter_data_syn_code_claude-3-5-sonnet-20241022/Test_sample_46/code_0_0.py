def can_execute_method(method, machines, parts):
    a, b, c = machines
    x, y, z = parts
    if method == 1:
        return a >= 1
    elif method == 2:
        return b >= 2
    elif method == 3:
        return c >= 2
    elif method == 4:
        return b >= 1 and c >= 1
    elif method == 5:
        return x >= 1 and y >= 1
    return False

def execute_method(method, machines, parts):
    a, b, c = machines
    x, y, z = parts
    
    if method == 1:  # A -> 2X + Y
        a -= 1
        x += 2
        y += 1
    elif method == 2:  # 2B -> X
        b -= 2
        x += 1
    elif method == 3:  # 2C -> Y
        c -= 2
        y += 1
    elif method == 4:  # B + C -> A
        b -= 1
        c -= 1
        a += 1
    elif method == 5:  # X + Y -> Z
        x -= 1
        y -= 1
        z += 1
    
    return [a, b, c], [x, y, z]

def simulate_dismantling():
    # Initial quantities
    machines = [2, 2, 2]  # A, B, C
    parts = [0, 0, 0]     # X, Y, Z
    
    current_method = 1
    changes_made = True
    
    while changes_made:
        changes_made = False
        for _ in range(5):  # Try each method once in a cycle
            if can_execute_method(current_method, machines, parts):
                machines, parts = execute_method(current_method, machines, parts)
                changes_made = True
            current_method = (current_method % 5) + 1
    
    result = machines + parts
    print([str(x) for x in result])

simulate_dismantling()