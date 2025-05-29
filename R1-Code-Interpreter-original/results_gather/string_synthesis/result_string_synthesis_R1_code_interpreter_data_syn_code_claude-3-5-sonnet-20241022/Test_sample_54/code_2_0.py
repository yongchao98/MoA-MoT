def apply_methods(blocks):
    new_states = []
    
    # Method 1: [A] + [B] + [C] -> {A}
    if blocks['[A]'] >= 1 and blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
        new_state = blocks.copy()
        new_state['[A]'] -= 1
        new_state['[B]'] -= 1
        new_state['[C]'] -= 1
        new_state['{A}'] += 1
        new_states.append(new_state)
    
    # Method 2: [A] + [B] -> {C}
    if blocks['[A]'] >= 1 and blocks['[B]'] >= 1:
        new_state = blocks.copy()
        new_state['[A]'] -= 1
        new_state['[B]'] -= 1
        new_state['{C}'] += 1
        new_states.append(new_state)
    
    # Method 3: [B] + [C] -> {B}
    if blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
        new_state = blocks.copy()
        new_state['[B]'] -= 1
        new_state['[C]'] -= 1
        new_state['{B}'] += 1
        new_states.append(new_state)
    
    # Method 4: [C] + [C] -> {C}
    if blocks['[C]'] >= 2:
        new_state = blocks.copy()
        new_state['[C]'] -= 2
        new_state['{C}'] += 1
        new_states.append(new_state)
    
    # Method 5: {A} + {C} -> (A) + (B)
    if blocks['{A}'] >= 1 and blocks['{C}'] >= 1:
        new_state = blocks.copy()
        new_state['{A}'] -= 1
        new_state['{C}'] -= 1
        new_state['(A)'] += 1
        new_state['(B)'] += 1
        new_states.append(new_state)
    
    # Method 6: {B} + {B} -> (C)
    if blocks['{B}'] >= 2:
        new_state = blocks.copy()
        new_state['{B}'] -= 2
        new_state['(C)'] += 1
        new_states.append(new_state)
    
    return new_states

def state_to_string(state):
    return ','.join(f"{k}:{v}" for k, v in sorted(state.items()) if v > 0)

def synthesize(blocks):
    current = blocks.copy()
    changed = True
    
    while changed:
        changed = False
        for _ in range(6):  # Try each method multiple times
            # Method 1
            while current['[A]'] >= 1 and current['[B]'] >= 1 and current['[C]'] >= 1:
                current['[A]'] -= 1
                current['[B]'] -= 1
                current['[C]'] -= 1
                current['{A}'] += 1
                changed = True
            
            # Method 2
            while current['[A]'] >= 1 and current['[B]'] >= 1:
                current['[A]'] -= 1
                current['[B]'] -= 1
                current['{C}'] += 1
                changed = True
            
            # Method 3
            while current['[B]'] >= 1 and current['[C]'] >= 1:
                current['[B]'] -= 1
                current['[C]'] -= 1
                current['{B}'] += 1
                changed = True
            
            # Method 4
            while current['[C]'] >= 2:
                current['[C]'] -= 2
                current['{C}'] += 1
                changed = True
            
            # Method 5
            while current['{A}'] >= 1 and current['{C}'] >= 1:
                current['{A}'] -= 1
                current['{C}'] -= 1
                current['(A)'] += 1
                current['(B)'] += 1
                changed = True
            
            # Method 6
            while current['{B}'] >= 2:
                current['{B}'] -= 2
                current['(C)'] += 1
                changed = True
    
    return current

initial_blocks = {
    '[A]': 2, '[B]': 7, '[C]': 3,
    '{A}': 0, '{B}': 0, '{C}': 0,
    '(A)': 0, '(B)': 0, '(C)': 0
}

result = synthesize(initial_blocks)
output = ', '.join(f"{k}: {v}" for k, v in sorted(result.items()) if v > 0)
print(output)