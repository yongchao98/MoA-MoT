def can_apply_method(method, machines, parts):
    A, B, C = machines
    X, Y, Z = parts
    if method == 1:
        return A >= 1
    elif method == 2:
        return B >= 2
    elif method == 3:
        return C >= 2
    elif method == 4:
        return B >= 1 and C >= 1
    elif method == 5:
        return X >= 1 and Y >= 1
    return False

def apply_method(method, machines, parts):
    A, B, C = machines
    X, Y, Z = parts
    
    if method == 1:
        A -= 1
        X += 2
        Y += 1
    elif method == 2:
        B -= 2
        X += 1
    elif method == 3:
        C -= 2
        Y += 1
    elif method == 4:
        B -= 1
        C -= 1
        A += 1
    elif method == 5:
        X -= 1
        Y -= 1
        Z += 1
    
    return [A, B, C], [X, Y, Z]

def simulate_dismantling():
    machines = [3, 2, 4]  # Initial A, B, C
    parts = [0, 0, 0]     # Initial X, Y, Z
    
    while True:
        any_method_applied = False
        
        for method in range(1, 6):
            while can_apply_method(method, machines, parts):
                machines, parts = apply_method(method, machines, parts)
                any_method_applied = True
        
        if not any_method_applied:
            break
    
    result = [str(x) for x in machines + parts]
    print(result)

simulate_dismantling()