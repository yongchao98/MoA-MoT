def apply_rule1(state):
    # [A] + [B] + [C] -> {A}
    if state['[A]'] >= 1 and state['[B]'] >= 1 and state['[C]'] >= 1:
        state['[A]'] -= 1
        state['[B]'] -= 1
        state['[C]'] -= 1
        state['{A}'] += 1
        return True
    return False

def apply_rule2(state):
    # [A] + [B] -> {C}
    if state['[A]'] >= 1 and state['[B]'] >= 1:
        state['[A]'] -= 1
        state['[B]'] -= 1
        state['{C}'] += 1
        return True
    return False

def apply_rule3(state):
    # [B] + [C] -> {B}
    if state['[B]'] >= 1 and state['[C]'] >= 1:
        state['[B]'] -= 1
        state['[C]'] -= 1
        state['{B}'] += 1
        return True
    return False

def apply_rule4(state):
    # [C] + [C] -> {C}
    if state['[C]'] >= 2:
        state['[C]'] -= 2
        state['{C}'] += 1
        return True
    return False

def apply_rule5(state):
    # {A} + {C} -> (A) + (B)
    if state['{A}'] >= 1 and state['{C}'] >= 1:
        state['{A}'] -= 1
        state['{C}'] -= 1
        state['(A)'] += 1
        state['(B)'] += 1
        return True
    return False

def apply_rule6(state):
    # {B} + {B} -> (C)
    if state['{B}'] >= 2:
        state['{B}'] -= 2
        state['(C)'] += 1
        return True
    return False

def synthesize():
    # Start from the current state
    state = {
        '[A]': 0, '[B]': 0, '[C]': 1,
        '{A}': 2, '{B}': 1, '{C}': 0,
        '(A)': 0, '(B)': 0, '(C)': 0
    }
    
    iterations = 0
    max_iterations = 100  # Safety limit
    
    while iterations < max_iterations:
        changes = False
        # Try each rule in order multiple times
        for rule in [apply_rule1, apply_rule2, apply_rule3, apply_rule4, apply_rule5, apply_rule6]:
            while rule(state):
                changes = True
        
        if not changes:
            break
        iterations += 1
    
    # Format the result
    result = []
    for block in ['[A]', '[B]', '[C]', '{A}', '{B}', '{C}', '(A)', '(B)', '(C)']:
        if state[block] > 0:
            result.append(f"{state[block]}{block}")
    
    print(" ".join(result))

synthesize()