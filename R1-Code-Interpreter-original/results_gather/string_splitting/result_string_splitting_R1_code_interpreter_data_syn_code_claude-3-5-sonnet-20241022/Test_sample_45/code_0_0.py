def can_execute_method(method, a, b, c, x, y, z):
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

def execute_method(method, a, b, c, x, y, z):
    if method == 1:  # A -> 2X + Y
        return a-1, b, c, x+2, y+1, z
    elif method == 2:  # 2B -> X
        return a, b-2, c, x+1, y, z
    elif method == 3:  # 2C -> Y
        return a, b, c-2, x, y+1, z
    elif method == 4:  # B + C -> A
        return a+1, b-1, c-1, x, y, z
    elif method == 5:  # X + Y -> Z
        return a, b, c, x-1, y-1, z+1
    return a, b, c, x, y, z

def solve_dismantling():
    # Initial quantities
    a, b, c = 4, 4, 3
    x, y, z = 0, 0, 0
    current_method = 1
    
    while True:
        # Try to find any executable method
        methods_tried = 0
        can_continue = False
        
        while methods_tried < 5:
            if can_execute_method(current_method, a, b, c, x, y, z):
                a, b, c, x, y, z = execute_method(current_method, a, b, c, x, y, z)
                can_continue = True
                break
            current_method = (current_method % 5) + 1
            methods_tried += 1
        
        if not can_continue:
            break
    
    return [str(a), str(b), str(c), str(x), str(y), str(z)]

result = solve_dismantling()
print(result)