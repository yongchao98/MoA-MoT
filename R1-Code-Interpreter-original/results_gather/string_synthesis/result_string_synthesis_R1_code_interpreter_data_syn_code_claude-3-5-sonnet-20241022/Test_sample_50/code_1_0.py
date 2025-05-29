def apply_rules(state):
    # Returns (rule_applied, new_state) if a rule can be applied, else (None, state)
    
    # Rule 1: [A] + [B] + [C] -> {A}
    if state['[A]'] >= 1 and state['[B]'] >= 1 and state['[C]'] >= 1:
        new_state = state.copy()
        new_state['[A]'] -= 1
        new_state['[B]'] -= 1
        new_state['[C]'] -= 1
        new_state['{A}'] += 1
        return (1, new_state)
    
    # Rule 2: [A] + [B] -> {C}
    if state['[A]'] >= 1 and state['[B]'] >= 1:
        new_state = state.copy()
        new_state['[A]'] -= 1
        new_state['[B]'] -= 1
        new_state['{C}'] += 1
        return (2, new_state)
    
    # Rule 3: [B] + [C] -> {B}
    if state['[B]'] >= 1 and state['[C]'] >= 1:
        new_state = state.copy()
        new_state['[B]'] -= 1
        new_state['[C]'] -= 1
        new_state['{B}'] += 1
        return (3, new_state)
    
    # Rule 4: [C] + [C] -> {C}
    if state['[C]'] >= 2:
        new_state = state.copy()
        new_state['[C]'] -= 2
        new_state['{C}'] += 1
        return (4, new_state)
    
    # Rule 5: {A} + {C} -> (A) + (B)
    if state['{A}'] >= 1 and state['{C}'] >= 1:
        new_state = state.copy()
        new_state['{A}'] -= 1
        new_state['{C}'] -= 1
        new_state['(A)'] += 1
        new_state['(B)'] += 1
        return (5, new_state)
    
    # Rule 6: {B} + {B} -> (C)
    if state['{B}'] >= 2:
        new_state = state.copy()
        new_state['{B}'] -= 2
        new_state['(C)'] += 1
        return (6, new_state)
    
    return (None, state)

def simulate():
    state = {
        '[A]': 4, '[B]': 5, '[C]': 5,
        '{A}': 0, '{B}': 0, '{C}': 0,
        '(A)': 0, '(B)': 0, '(C)': 0
    }
    
    step = 1
    while True:
        rule, new_state = apply_rules(state)
        if rule is None:
            break
        state = new_state
        print(f"Step {step}, Rule {rule}:", end=" ")
        result = ""
        for block in ['[A]', '[B]', '[C]', '{A}', '{B}', '{C}', '(A)', '(B)', '(C)']:
            if state[block] > 0:
                result += f"{state[block]} {block}, "
        print(result[:-2])
        step += 1
    
    print("\nFinal state:", end=" ")
    result = ""
    for block in ['[A]', '[B]', '[C]', '{A}', '{B}', '{C}', '(A)', '(B)', '(C)']:
        if state[block] > 0:
            result += f"{state[block]} {block}, "
    print(result[:-2])

simulate()