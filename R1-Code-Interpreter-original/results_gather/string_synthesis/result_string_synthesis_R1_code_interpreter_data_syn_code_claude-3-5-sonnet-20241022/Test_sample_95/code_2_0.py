def apply_rules_until_stable(state):
    while True:
        old_state = state.copy()
        
        # Rule 4: [C] + [C] -> {C}
        pairs = state['[C]'] // 2
        if pairs > 0:
            state['[C]'] -= 2 * pairs
            state['{C}'] += pairs
            
        # Rule 2: [A] + [B] -> {C}
        pairs = min(state['[A]'], state['[B]'])
        if pairs > 0:
            state['[A]'] -= pairs
            state['[B]'] -= pairs
            state['{C}'] += pairs
            
        # Rule 1: [A] + [B] + [C] -> {A}
        triples = min(state['[A]'], state['[B]'], state['[C]'])
        if triples > 0:
            state['[A]'] -= triples
            state['[B]'] -= triples
            state['[C]'] -= triples
            state['{A}'] += triples
            
        # Rule 3: [B] + [C] -> {B}
        pairs = min(state['[B]'], state['[C]'])
        if pairs > 0:
            state['[B]'] -= pairs
            state['[C]'] -= pairs
            state['{B}'] += pairs
            
        # Rule 5: {A} + {C} -> (A) + (B)
        pairs = min(state['{A}'], state['{C}'])
        if pairs > 0:
            state['{A}'] -= pairs
            state['{C}'] -= pairs
            state['(A)'] += pairs
            state['(B)'] += pairs
            
        # Rule 6: {B} + {B} -> (C)
        pairs = state['{B}'] // 2
        if pairs > 0:
            state['{B}'] -= 2 * pairs
            state['(C)'] += pairs
            
        if state == old_state:
            break

state = {
    '[A]': 4, '[B]': 7, '[C]': 6,
    '{A}': 0, '{B}': 0, '{C}': 0,
    '(A)': 0, '(B)': 0, '(C)': 0
}

apply_rules_until_stable(state)

result = []
for block in ['[A]', '[B]', '[C]', '{A}', '{B}', '{C}', '(A)', '(B)', '(C)']:
    if state[block] > 0:
        result.append(f"{state[block]} {block}")
print(', '.join(result))