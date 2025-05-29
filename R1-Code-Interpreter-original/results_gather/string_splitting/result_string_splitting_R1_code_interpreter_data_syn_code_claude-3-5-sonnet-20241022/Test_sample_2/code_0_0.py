def can_apply_method(method, machines, parts):
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

def apply_method(method, machines, parts):
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

def solve_dismantling():
    # Initial state
    machines = [2, 2, 1]  # A, B, C
    parts = [0, 0, 0]    # X, Y, Z
    
    current_method = 1
    changes = True
    
    while changes:
        changes = False
        for _ in range(5):  # Try each method once per cycle
            if can_apply_method(current_method, machines, parts):
                machines, parts = apply_method(current_method, machines, parts)
                changes = True
            current_method = (current_method % 5) + 1
    
    # Convert to strings for output format
    result = [str(x) for x in machines + parts]
    print(result)

solve_dismantling()