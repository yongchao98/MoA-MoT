def apply_rules(state):
    # Keep track of initial state
    prev_state = state.copy()
    
    # Rule 4: [C] + [C] -> {C}
    while state['[C]'] >= 2:
        state['[C]'] -= 2
        state['{C}'] += 1
    
    # Rule 3: [B] + [C] -> {B}
    while state['[B]'] >= 1 and state['[C]'] >= 1:
        state['[B]'] -= 1
        state['[C]'] -= 1
        state['{B}'] += 1
    
    # Rule 2: [A] + [B] -> {C}
    while state['[A]'] >= 1 and state['[B]'] >= 1:
        state['[A]'] -= 1
        state['[B]'] -= 1
        state['{C}'] += 1
    
    # Rule 5: {A} + {C} -> (A) + (B)
    while state['{A}'] >= 1 and state['{C}'] >= 1:
        state['{A}'] -= 1
        state['{C}'] -= 1
        state['(A)'] += 1
        state['(B)'] += 1
    
    # Rule 6: {B} + {B} -> (C)
    while state['{B}'] >= 2:
        state['{B}'] -= 2
        state['(C)'] += 1
    
    return state != prev_state

def simulate():
    state = {
        '[A]': 5, '[B]': 2, '[C]': 7,
        '{A}': 0, '{B}': 0, '{C}': 0,
        '(A)': 0, '(B)': 0, '(C)': 0
    }
    
    while apply_rules(state):
        pass
    
    result = []
    for block in ['[A]', '[B]', '[C]', '{A}', '{B}', '{C}', '(A)', '(B)', '(C)']:
        if state[block] > 0:
            result.append(f"{state[block]} {block}")
    
    print(", ".join(result))

simulate()