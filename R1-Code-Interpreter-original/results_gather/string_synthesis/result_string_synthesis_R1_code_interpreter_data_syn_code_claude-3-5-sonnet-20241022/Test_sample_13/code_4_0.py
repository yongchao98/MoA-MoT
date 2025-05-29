def apply_rules(state):
    changed = True
    while changed:
        changed = False
        initial = state.copy()
        
        # Rule 4 first (fundamental transformation)
        if state['[C]'] >= 2:
            state['[C]'] -= 2
            state['{C}'] += 1
            changed = True
            continue
            
        # Rule 1 (creates {A} which is needed for Rule 5)
        if state['[A]'] >= 1 and state['[B]'] >= 1 and state['[C]'] >= 1:
            state['[A]'] -= 1
            state['[B]'] -= 1
            state['[C]'] -= 1
            state['{A}'] += 1
            changed = True
            continue
            
        # Rule 2 (creates {C})
        if state['[A]'] >= 1 and state['[B]'] >= 1:
            state['[A]'] -= 1
            state['[B]'] -= 1
            state['{C}'] += 1
            changed = True
            continue
            
        # Rule 3 (creates {B})
        if state['[B]'] >= 1 and state['[C]'] >= 1:
            state['[B]'] -= 1
            state['[C]'] -= 1
            state['{B}'] += 1
            changed = True
            continue
            
        # Rule 5 (creates highest order blocks)
        if state['{A}'] >= 1 and state['{C}'] >= 1:
            state['{A}'] -= 1
            state['{C}'] -= 1
            state['(A)'] += 1
            state['(B)'] += 1
            changed = True
            continue
            
        # Rule 6
        if state['{B}'] >= 2:
            state['{B}'] -= 2
            state['(C)'] += 1
            changed = True
            continue

def simulate():
    state = {
        '[A]': 5, '[B]': 2, '[C]': 7,
        '{A}': 0, '{B}': 0, '{C}': 0,
        '(A)': 0, '(B)': 0, '(C)': 0
    }
    
    apply_rules(state)
    
    result = []
    for block in ['[A]', '[B]', '[C]', '{A}', '{B}', '{C}', '(A)', '(B)', '(C)']:
        if state[block] > 0:
            result.append(f"{state[block]} {block}")
    
    print(", ".join(result))

simulate()