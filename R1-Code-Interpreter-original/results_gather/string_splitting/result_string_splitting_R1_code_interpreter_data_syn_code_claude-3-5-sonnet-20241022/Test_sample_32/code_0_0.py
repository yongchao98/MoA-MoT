def method1(inventory):
    # A -> 2X + Y
    if inventory['A'] >= 1:
        inventory['A'] -= 1
        inventory['X'] += 2
        inventory['Y'] += 1
        return True
    return False

def method2(inventory):
    # 2B -> X
    if inventory['B'] >= 2:
        inventory['B'] -= 2
        inventory['X'] += 1
        return True
    return False

def method3(inventory):
    # 2C -> Y
    if inventory['C'] >= 2:
        inventory['C'] -= 2
        inventory['Y'] += 1
        return True
    return False

def method4(inventory):
    # B + C -> A
    if inventory['B'] >= 1 and inventory['C'] >= 1:
        inventory['B'] -= 1
        inventory['C'] -= 1
        inventory['A'] += 1
        return True
    return False

def method5(inventory):
    # X + Y -> Z
    if inventory['X'] >= 1 and inventory['Y'] >= 1:
        inventory['X'] -= 1
        inventory['Y'] -= 1
        inventory['Z'] += 1
        return True
    return False

# Initialize inventory
inventory = {
    'A': 3,
    'B': 2,
    'C': 3,
    'X': 0,
    'Y': 0,
    'Z': 0
}

methods = [method1, method2, method3, method4, method5]
current_method = 0
changes_made = True

while changes_made:
    changes_made = False
    for i in range(len(methods)):
        if methods[current_method](inventory):
            changes_made = True
        current_method = (current_method + 1) % 5
    if not changes_made:
        break

# Format the result
result = [str(inventory['A']), str(inventory['B']), str(inventory['C']), 
          str(inventory['X']), str(inventory['Y']), str(inventory['Z'])]
print(result)