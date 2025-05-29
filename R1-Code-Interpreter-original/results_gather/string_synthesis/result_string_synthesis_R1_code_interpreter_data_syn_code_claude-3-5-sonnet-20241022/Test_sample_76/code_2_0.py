from collections import defaultdict

def state_to_string(state):
    return ' '.join(f"{v}{k}" for k, v in sorted(state.items()) if v > 0)

def apply_rules(state):
    states_seen = set()
    states_to_process = [state.copy()]
    final_states = set()
    
    while states_to_process:
        current_state = states_to_process.pop()
        state_str = state_to_string(current_state)
        
        if state_str in states_seen:
            continue
            
        states_seen.add(state_str)
        changes_made = False
        
        # Try Rule 1: [A] + [B] + [C] -> {A}
        if current_state['[A]'] >= 1 and current_state['[B]'] >= 1 and current_state['[C]'] >= 1:
            new_state = current_state.copy()
            new_state['[A]'] -= 1
            new_state['[B]'] -= 1
            new_state['[C]'] -= 1
            new_state['{A}'] += 1
            states_to_process.append(new_state)
            changes_made = True
            
        # Try Rule 2: [A] + [B] -> {C}
        if current_state['[A]'] >= 1 and current_state['[B]'] >= 1:
            new_state = current_state.copy()
            new_state['[A]'] -= 1
            new_state['[B]'] -= 1
            new_state['{C}'] += 1
            states_to_process.append(new_state)
            changes_made = True
            
        # Try Rule 3: [B] + [C] -> {B}
        if current_state['[B]'] >= 1 and current_state['[C]'] >= 1:
            new_state = current_state.copy()
            new_state['[B]'] -= 1
            new_state['[C]'] -= 1
            new_state['{B}'] += 1
            states_to_process.append(new_state)
            changes_made = True
            
        # Try Rule 4: [C] + [C] -> {C}
        if current_state['[C]'] >= 2:
            new_state = current_state.copy()
            new_state['[C]'] -= 2
            new_state['{C}'] += 1
            states_to_process.append(new_state)
            changes_made = True
            
        # Try Rule 5: {A} + {C} -> (A) + (B)
        if current_state['{A}'] >= 1 and current_state['{C}'] >= 1:
            new_state = current_state.copy()
            new_state['{A}'] -= 1
            new_state['{C}'] -= 1
            new_state['(A)'] += 1
            new_state['(B)'] += 1
            states_to_process.append(new_state)
            changes_made = True
            
        # Try Rule 6: {B} + {B} -> (C)
        if current_state['{B}'] >= 2:
            new_state = current_state.copy()
            new_state['{B}'] -= 2
            new_state['(C)'] += 1
            states_to_process.append(new_state)
            changes_made = True
            
        if not changes_made:
            final_states.add(state_str)
    
    return final_states

initial_state = {
    '[A]': 0, '[B]': 0, '[C]': 1,
    '{A}': 2, '{B}': 1, '{C}': 0,
    '(A)': 0, '(B)': 0, '(C)': 0
}

final_states = apply_rules(initial_state)
print("Final possible states:")
for state in sorted(final_states):
    print(state)