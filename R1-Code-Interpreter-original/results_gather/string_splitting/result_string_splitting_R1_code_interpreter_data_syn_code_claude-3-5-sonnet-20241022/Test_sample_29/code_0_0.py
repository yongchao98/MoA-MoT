def method1(state):
    # A -> 2X + Y
    if state['A'] >= 1:
        state['A'] -= 1
        state['X'] += 2
        state['Y'] += 1
        return True
    return False

def method2(state):
    # 2B -> X
    if state['B'] >= 2:
        state['B'] -= 2
        state['X'] += 1
        return True
    return False

def method3(state):
    # 2C -> Y
    if state['C'] >= 2:
        state['C'] -= 2
        state['Y'] += 1
        return True
    return False

def method4(state):
    # B + C -> A
    if state['B'] >= 1 and state['C'] >= 1:
        state['B'] -= 1
        state['C'] -= 1
        state['A'] += 1
        return True
    return False

def method5(state):
    # X + Y -> Z
    if state['X'] >= 1 and state['Y'] >= 1:
        state['X'] -= 1
        state['Y'] -= 1
        state['Z'] += 1
        return True
    return False

# Initial state
state = {
    'A': 2,
    'B': 1,
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
    for i in range(5):
        method_index = (current_method + i) % 5
        if methods[method_index](state):
            current_method = (method_index + 1) % 5
            changes_made = True
            break

# Format the result
result = [str(state['A']), str(state['B']), str(state['C']), 
          str(state['X']), str(state['Y']), str(state['Z'])]
print(result)