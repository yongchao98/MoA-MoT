def can_execute_method(inventory, method_num):
    A, B, C, X, Y, Z = inventory
    if method_num == 1:
        return A >= 1
    elif method_num == 2:
        return B >= 2
    elif method_num == 3:
        return C >= 2
    elif method_num == 4:
        return B >= 1 and C >= 1
    elif method_num == 5:
        return X >= 1 and Y >= 1
    return False

def execute_method(inventory, method_num):
    A, B, C, X, Y, Z = inventory
    if method_num == 1:  # A -> 2X + Y
        A -= 1
        X += 2
        Y += 1
    elif method_num == 2:  # 2B -> X
        B -= 2
        X += 1
    elif method_num == 3:  # 2C -> Y
        C -= 2
        Y += 1
    elif method_num == 4:  # B + C -> A
        B -= 1
        C -= 1
        A += 1
    elif method_num == 5:  # X + Y -> Z
        X -= 1
        Y -= 1
        Z += 1
    return [A, B, C, X, Y, Z]

# Initial inventory: [A, B, C, X, Y, Z]
inventory = [2, 3, 2, 0, 0, 0]
current_method = 1
changes_made = True

print("Initial state:", inventory)

while changes_made:
    changes_made = False
    for _ in range(5):  # Try each method once in a cycle
        if can_execute_method(inventory, current_method):
            inventory = execute_method(inventory, current_method)
            changes_made = True
            print(f"After method {current_method}:", inventory)
        current_method = (current_method % 5) + 1

print("\nFinal result:", inventory)