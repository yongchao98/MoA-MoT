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
    # Initial quantities
    machines = [2, 3, 3]  # A, B, C
    parts = [0, 0, 0]    # X, Y, Z
    
    current_method = 1
    no_method_executable = False
    
    while not no_method_executable:
        # Try one complete cycle
        methods_tried = 0
        method_executed = False
        
        while methods_tried < 5:
            if can_execute_method(current_method, machines, parts):
                machines, parts = execute_method(current_method, machines, parts)
                method_executed = True
                break
            
            current_method = (current_method % 5) + 1
            methods_tried += 1
        
        if not method_executed:
            no_method_executable = True
    
    result = machines + parts
    print([str(x) for x in result])

solve_dismantling()