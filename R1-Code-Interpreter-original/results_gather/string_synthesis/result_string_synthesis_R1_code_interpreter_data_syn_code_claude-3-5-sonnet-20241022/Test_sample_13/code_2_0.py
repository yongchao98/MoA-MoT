def apply_rules(state):
    changes_made = False
    
    # Rule 1: [A] + [B] + [C] -> {A}
    if state['[A]'] >= 1 and state['[B]'] >= 1 and state['[C]'] >= 1:
        state['[A]'] -= 1
        state['[B]'] -= 1
        state['[C]'] -= 1
        state['{A}'] += 1
        changes_made = True
        return True

    # Rule 2: [A] + [B] -> {C}
    if state['[A]'] >= 1 and state['[B]'] >= 1:
        state['[A]'] -= 1
        state['[B]'] -= 1
        state['{C}'] += 1
        changes_made = True
        return True

    # Rule 3: [B] + [C] -> {B}
    if state['[B]'] >= 1 and state['[C]'] >= 1:
        state['[B]'] -= 1
        state['[C]'] -= 1
        state['{B}'] += 1
        changes_made = True
        return True

    # Rule 4: [C] + [C] -> {C}
    if state['[C]'] >= 2:
        state['[C]'] -= 2
        state['{C}'] += 1
        changes_made = True
        return True

    # Rule 5: {A} + {C} -> (A) + (B)
    if state['{A}'] >= 1 and state['{C}'] >= 1:
        state['{A}'] -= 1
        state['{C}'] -= 1
        state['(A)'] += 1
        state['(B)'] += 1
        changes_made = True
        return True

    # Rule 6: {B} + {B} -> (C)
    if state['{B}'] >= 2:
        state['{B}'] -= 2
        state['(C)'] += 1
        changes_made = True
        return True

    return False

def simulate():
    state = {
        '[A]': 5, '[B]': 2, '[C]': 7,
        '{A}': 0, '{B}': 0, '{C}': 0,
        '(A)': 0, '(B)': 0, '(C)': 0
    }
    
    # Keep applying rules until no more changes
    while True:
        if not apply_rules(state):
            break
    
    # Print final state
    result = []
    for block in ['[A]', '[B]', '[C]', '{A}', '{B}', '{C}', '(A)', '(B)', '(C)']:
        if state[block] > 0:
            result.append(f"{state[block]} {block}")
    
    print(", ".join(result))

simulate()