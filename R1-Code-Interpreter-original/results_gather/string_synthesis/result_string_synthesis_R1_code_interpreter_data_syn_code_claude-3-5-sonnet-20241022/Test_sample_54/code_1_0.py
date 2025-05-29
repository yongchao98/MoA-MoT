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

def synthesize(initial_blocks):
    seen_states = set()
    states_to_check = [initial_blocks]
    final_states = set()
    
    while states_to_check:
        current_state = states_to_check.pop(0)
        state_str = state_to_string(current_state)
        
        if state_str in seen_states:
            continue
            
        seen_states.add(state_str)
        new_states = apply_methods(current_state)
        
        if not new_states:
            final_states.add(state_str)
        else:
            states_to_check.extend(new_states)
    
    # Find the state with the most higher-tier blocks
    best_state = max(final_states, key=lambda x: sum(int(v) for k,v in dict(item.split(':') for item in x.split(',')).items() if '{' in k or '(' in k))
    
    return dict(item.split(':') for item in best_state.split(','))

initial_blocks = {
    '[A]': 2, '[B]': 7, '[C]': 3,
    '{A}': 0, '{B}': 0, '{C}': 0,
    '(A)': 0, '(B)': 0, '(C)': 0
}

result = synthesize(initial_blocks)
print(', '.join(f"{k}: {v}" for k, v in sorted(result.items()) if v > 0))