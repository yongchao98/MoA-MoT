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
    return a, b, c, x, y, z

def solve_dismantling():
    # Initial quantities
    a, b, c = 4, 4, 4
    x, y, z = 0, 0, 0
    
    current_method = 1
    made_change = True
    
    while made_change:
        made_change = False
        for _ in range(5):  # Try each method once in a cycle
            if can_execute_method(current_method, a, b, c, x, y, z):
                a, b, c, x, y, z = execute_method(current_method, a, b, c, x, y, z)
                made_change = True
            current_method = (current_method % 5) + 1
    
    return [str(a), str(b), str(c), str(x), str(y), str(z)]

result = solve_dismantling()
print(result)