def can_execute_method(method, inventory):
    A, B, C, X, Y, Z = inventory
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

def execute_method(method, inventory):
    A, B, C, X, Y, Z = inventory
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
    return [A, B, C, X, Y, Z]

def simulate_dismantling():
    # Initial inventory [A, B, C, X, Y, Z]
    inventory = [3, 2, 2, 0, 0, 0]
    current_method = 1
    changes_made = True
    
    while changes_made:
        changes_made = False
        for _ in range(5):  # Try each method once in a cycle
            if can_execute_method(current_method, inventory):
                inventory = execute_method(current_method, inventory)
                changes_made = True
            current_method = (current_method % 5) + 1
    
    return [str(x) for x in inventory]

# Run simulation and print result
print(simulate_dismantling())