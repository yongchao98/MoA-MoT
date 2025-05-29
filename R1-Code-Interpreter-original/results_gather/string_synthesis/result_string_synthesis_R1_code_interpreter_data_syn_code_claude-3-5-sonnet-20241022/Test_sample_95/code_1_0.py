def apply_rules(state):
    changes = True
    while changes:
        changes = False
        # Rule 4: [C] + [C] -> {C}
        while state['[C]'] >= 2:
            state['[C]'] -= 2
            state['{C}'] += 1
            changes = True
        
        # Rule 2: [A] + [B] -> {C}
        while state['[A]'] >= 1 and state['[B]'] >= 1:
            state['[A]'] -= 1
            state['[B]'] -= 1
            state['{C}'] += 1
            changes = True
            
        # Rule 1: [A] + [B] + [C] -> {A}
        while state['[A]'] >= 1 and state['[B]'] >= 1 and state['[C]'] >= 1:
            state['[A]'] -= 1
            state['[B]'] -= 1
            state['[C]'] -= 1
            state['{A}'] += 1
            changes = True
            
        # Rule 3: [B] + [C] -> {B}
        while state['[B]'] >= 1 and state['[C]'] >= 1:
            state['[B]'] -= 1
            state['[C]'] -= 1
            state['{B}'] += 1
            changes = True
            
        # Rule 5: {A} + {C} -> (A) + (B)
        while state['{A}'] >= 1 and state['{C}'] >= 1:
            state['{A}'] -= 1
            state['{C}'] -= 1
            state['(A)'] += 1
            state['(B)'] += 1
            changes = True
            
        # Rule 6: {B} + {B} -> (C)
        while state['{B}'] >= 2:
            state['{B}'] -= 2
            state['(C)'] += 1
            changes = True

state = {
    '[A]': 4, '[B]': 7, '[C]': 6,
    '{A}': 0, '{B}': 0, '{C}': 0,
    '(A)': 0, '(B)': 0, '(C)': 0
}

apply_rules(state)

result = []
for block in ['[A]', '[B]', '[C]', '{A}', '{B}', '{C}', '(A)', '(B)', '(C)']:
    if state[block] > 0:
        result.append(f"{state[block]} {block}")
print(', '.join(result))