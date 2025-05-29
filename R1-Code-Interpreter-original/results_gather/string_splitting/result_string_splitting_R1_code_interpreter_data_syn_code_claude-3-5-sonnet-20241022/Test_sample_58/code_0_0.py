def can_execute_method(method, machines, parts):
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

def execute_method(method, machines, parts):
    A, B, C = machines
    X, Y, Z = parts
    
    if method == 1:  # A -> 2X + Y
        A -= 1
        X += 2
        Y += 1
    elif method == 2:  # 2B -> X
        B -= 2
        X += 1
    elif method == 3:  # 2C -> Y
        C -= 2
        Y += 1
    elif method == 4:  # B + C -> A
        B -= 1
        C -= 1
        A += 1
    elif method == 5:  # X + Y -> Z
        X -= 1
        Y -= 1
        Z += 1
    
    return [A, B, C], [X, Y, Z]

def simulate_dismantling():
    # Initial quantities
    machines = [4, 2, 4]  # A, B, C
    parts = [0, 0, 0]    # X, Y, Z
    
    current_method = 1
    can_continue = True
    
    while can_continue:
        can_continue = False
        # Try all methods in sequence
        for _ in range(5):
            if can_execute_method(current_method, machines, parts):
                machines, parts = execute_method(current_method, machines, parts)
                can_continue = True
            
            current_method = (current_method % 5) + 1
    
    result = [str(x) for x in machines + parts]
    print(result)

simulate_dismantling()