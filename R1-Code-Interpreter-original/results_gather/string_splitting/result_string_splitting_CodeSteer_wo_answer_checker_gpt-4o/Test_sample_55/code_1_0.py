# Initialize the quantities
machines = {'A': 2, 'B': 4, 'C': 2}
parts = {'X': 0, 'Y': 0, 'Z': 0}

# Define the methods
def method1():
    if machines['A'] >= 1:
        machines['A'] -= 1
        parts['X'] += 2
        parts['Y'] += 1
        return True
    return False

def method2():
    if machines['B'] >= 2:
        machines['B'] -= 2
        parts['X'] += 1
        return True
    return False

def method3():
    if machines['C'] >= 2:
        machines['C'] -= 2
        parts['Y'] += 1
        return True
    return False

def method4():
    if machines['B'] >= 1 and machines['C'] >= 1:
        machines['B'] -= 1
        machines['C'] -= 1
        machines['A'] += 1
        return True
    return False

def method5():
    if parts['X'] >= 1 and parts['Y'] >= 1:
        parts['X'] -= 1
        parts['Y'] -= 1
        parts['Z'] += 1
        return True
    return False

# Cycle through the methods
methods = [method1, method2, method3, method4, method5]
while True:
    changes = False
    for method in methods:
        if method():
            changes = True
    if not changes:
        break

# Output the final quantities
result = [str(machines['A']), str(machines['B']), str(machines['C']),
          str(parts['X']), str(parts['Y']), str(parts['Z'])]
print(result)