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
    'A': 1,
    'B': 3,
    'C': 3,
    'X': 0,
    'Y': 0,
    'Z': 0
}

methods = [method1, method2, method3, method4, method5]
method_index = 0
changes_made = True

while changes_made:
    changes_made = False
    for _ in range(len(methods)):
        if methods[method_index](state):
            changes_made = True
        method_index = (method_index + 1) % len(methods)
    if not changes_made:
        break

result = [str(state['A']), str(state['B']), str(state['C']), 
          str(state['X']), str(state['Y']), str(state['Z'])]
print(result)