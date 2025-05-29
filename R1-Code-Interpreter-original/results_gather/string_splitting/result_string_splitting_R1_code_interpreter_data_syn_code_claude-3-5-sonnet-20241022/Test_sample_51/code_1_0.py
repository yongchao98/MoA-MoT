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

def solve_dismantling():
    machines = [2, 3, 4]  # A, B, C
    parts = [0, 0, 0]    # X, Y, Z
    current_method = 1
    
    while True:
        if can_execute_method(current_method, machines, parts):
            machines, parts = execute_method(current_method, machines, parts)
        current_method = current_method % 5 + 1
        
        # Check if any method can be executed
        can_continue = False
        for m in range(1, 6):
            if can_execute_method(m, machines, parts):
                can_continue = True
                break
        if not can_continue:
            break
    
    print([str(x) for x in machines + parts])

solve_dismantling()