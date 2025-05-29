def can_execute_method(method, quantities):
    a, b, c, x, y, z = quantities
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

def execute_method(method, quantities):
    a, b, c, x, y, z = quantities
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

def simulate_dismantling():
    # Initial quantities [A, B, C, X, Y, Z]
    quantities = [1, 1, 3, 0, 0, 0]
    current_method = 1
    no_method_possible = False
    
    while not no_method_possible:
        no_method_possible = True
        # Try each method in sequence
        for i in range(5):
            method = ((current_method + i - 1) % 5) + 1
            if can_execute_method(method, quantities):
                quantities = execute_method(method, quantities)
                current_method = (method % 5) + 1
                no_method_possible = False
                break
    
    return [str(q) for q in quantities]

result = simulate_dismantling()
print(result)