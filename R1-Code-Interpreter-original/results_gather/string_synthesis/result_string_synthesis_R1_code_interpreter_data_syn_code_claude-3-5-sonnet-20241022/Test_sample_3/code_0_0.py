def apply_rule1(state):
    # [A] + [B] + [C] -> {A}
    if state['[A]'] >= 1 and state['[B]'] >= 1 and state['[C]'] >= 1:
        new_state = state.copy()
        new_state['[A]'] -= 1
        new_state['[B]'] -= 1
        new_state['[C]'] -= 1
        new_state['{A}'] += 1
        return new_state
    return None

def apply_rule2(state):
    # [A] + [B] -> {C}
    if state['[A]'] >= 1 and state['[B]'] >= 1:
        new_state = state.copy()
        new_state['[A]'] -= 1
        new_state['[B]'] -= 1
        new_state['{C}'] += 1
        return new_state
    return None

def apply_rule3(state):
    # [B] + [C] -> {B}
    if state['[B]'] >= 1 and state['[C]'] >= 1:
        new_state = state.copy()
        new_state['[B]'] -= 1
        new_state['[C]'] -= 1
        new_state['{B}'] += 1
        return new_state
    return None

def apply_rule4(state):
    # [C] + [C] -> {C}
    if state['[C]'] >= 2:
        new_state = state.copy()
        new_state['[C]'] -= 2
        new_state['{C}'] += 1
        return new_state
    return None

def apply_rule5(state):
    # {A} + {C} -> (A) + (B)
    if state['{A}'] >= 1 and state['{C}'] >= 1:
        new_state = state.copy()
        new_state['{A}'] -= 1
        new_state['{C}'] -= 1
        new_state['(A)'] += 1
        new_state['(B)'] += 1
        return new_state
    return None

def apply_rule6(state):
    # {B} + {B} -> (C)
    if state['{B}'] >= 2:
        new_state = state.copy()
        new_state['{B}'] -= 2
        new_state['(C)'] += 1
        return new_state
    return None

def process_synthesis():
    # Initial state
    state = {
        '[A]': 5, '[B]': 1, '[C]': 1,
        '{A}': 0, '{B}': 0, '{C}': 0,
        '(A)': 0, '(B)': 0, '(C)': 0
    }
    
    rules = [apply_rule1, apply_rule2, apply_rule3, apply_rule4, apply_rule5, apply_rule6]
    changes = True
    
    while changes:
        changes = False
        for rule in rules:
            new_state = rule(state)
            if new_state:
                state = new_state
                changes = True
                break
    
    # Format result
    result = []
    for block, count in state.items():
        if count > 0:
            result.append(f"{count} {block}")
    
    print(", ".join(result))

process_synthesis()