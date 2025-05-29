def method1(resources):
    # A -> 2X + Y
    if resources['A'] >= 1:
        resources['A'] -= 1
        resources['X'] += 2
        resources['Y'] += 1
        return True
    return False

def method2(resources):
    # 2B -> X
    if resources['B'] >= 2:
        resources['B'] -= 2
        resources['X'] += 1
        return True
    return False

def method3(resources):
    # 2C -> Y
    if resources['C'] >= 2:
        resources['C'] -= 2
        resources['Y'] += 1
        return True
    return False

def method4(resources):
    # B + C -> A
    if resources['B'] >= 1 and resources['C'] >= 1:
        resources['B'] -= 1
        resources['C'] -= 1
        resources['A'] += 1
        return True
    return False

def method5(resources):
    # X + Y -> Z
    if resources['X'] >= 1 and resources['Y'] >= 1:
        resources['X'] -= 1
        resources['Y'] -= 1
        resources['Z'] += 1
        return True
    return False

# Initialize resources
resources = {
    'A': 3,
    'B': 3,
    'C': 2,
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
        if methods[current_method](resources):
            changes_made = True
        current_method = (current_method + 1) % len(methods)
    if not changes_made:
        break

result = [str(resources['A']), str(resources['B']), str(resources['C']), 
          str(resources['X']), str(resources['Y']), str(resources['Z'])]
print(result)