# Initial quantities
machines = {'A': 2, 'B': 4, 'C': 4}
parts = {'X': 0, 'Y': 0, 'Z': 0}

# Define the methods
def method_1():
    if machines['A'] >= 1:
        machines['A'] -= 1
        parts['X'] += 2
        parts['Y'] += 1
        return True
    return False

def method_2():
    if machines['B'] >= 2:
        machines['B'] -= 2
        parts['X'] += 1
        return True
    return False

def method_3():
    if machines['C'] >= 2:
        machines['C'] -= 2
        parts['Y'] += 1
        return True
    return False

def method_4():
    if machines['B'] >= 1 and machines['C'] >= 1:
        machines['B'] -= 1
        machines['C'] -= 1
        machines['A'] += 1
        return True
    return False

def method_5():
    if parts['X'] >= 1 and parts['Y'] >= 1:
        parts['X'] -= 1
        parts['Y'] -= 1
        parts['Z'] += 1
        return True
    return False

# Cycle through the methods
methods = [method_1, method_2, method_3, method_4, method_5]
while True:
    executed = False
    for method in methods:
        if method():
            executed = True
            break
    if not executed:
        break

# Output the final quantities
result = [str(machines['A']), str(machines['B']), str(machines['C']), str(parts['X']), str(parts['Y']), str(parts['Z'])]
print(result)