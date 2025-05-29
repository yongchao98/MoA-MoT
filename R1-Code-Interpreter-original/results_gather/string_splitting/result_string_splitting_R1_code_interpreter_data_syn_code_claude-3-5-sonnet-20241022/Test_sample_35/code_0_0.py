def can_execute_method(method, inventory):
    a, b, c, x, y, z = inventory
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

def execute_method(method, inventory):
    a, b, c, x, y, z = inventory
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
    return [a, b, c, x, y, z]

def solve_dismantling():
    # Initial inventory [A, B, C, X, Y, Z]
    inventory = [3, 3, 3, 0, 0, 0]
    current_method = 1
    
    while True:
        # Try to find any executable method
        methods_tried = 0
        while methods_tried < 5:
            if can_execute_method(current_method, inventory):
                inventory = execute_method(current_method, inventory)
                break
            current_method = (current_method % 5) + 1
            methods_tried += 1
        
        # If no method can be executed, we're done
        if methods_tried == 5:
            break
        
        # Move to next method
        current_method = (current_method % 5) + 1
    
    # Convert results to strings
    return [str(x) for x in inventory]

# Run the simulation and print the result
print(solve_dismantling())