def can_execute(method, machines, parts):
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
    
    if method == 1:
        a -= 1
        x += 2
        y += 1
    elif method == 2:
        b -= 2
        x += 1
    elif method == 3:
        c -= 2
        y += 1
    elif method == 4:
        b -= 1
        c -= 1
        a += 1
    elif method == 5:
        x -= 1
        y -= 1
        z += 1
    
    return [a, b, c], [x, y, z]

def solve_dismantling():
    machines = [2, 3, 2]  # A, B, C
    parts = [0, 0, 0]     # X, Y, Z
    current_method = 1
    made_progress = True
    
    while made_progress:
        made_progress = False
        for _ in range(5):  # Try each method once in a cycle
            if can_execute(current_method, machines, parts):
                machines, parts = execute_method(current_method, machines, parts)
                made_progress = True
            current_method = (current_method % 5) + 1
    
    return machines + parts

# Run simulation and format output
result = solve_dismantling()
print([str(x) for x in result])