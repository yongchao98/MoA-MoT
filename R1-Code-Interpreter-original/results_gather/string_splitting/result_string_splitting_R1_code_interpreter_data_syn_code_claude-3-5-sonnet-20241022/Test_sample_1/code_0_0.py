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

def simulate():
    # Initial quantities
    machines = [3, 1, 3]  # A, B, C
    parts = [0, 0, 0]    # X, Y, Z
    
    while True:
        progress = False
        # Try each method in order
        for method in range(1, 6):
            while can_execute_method(method, machines, parts):
                machines, parts = execute_method(method, machines, parts)
                progress = True
        
        if not progress:
            break
    
    return machines + parts

result = simulate()
print([str(x) for x in result])