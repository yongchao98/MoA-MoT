def print_state(state, step=""):
    blocks = [k for k,v in state.items() if v > 0]
    result = [f"{state[b]} {b}" for b in blocks]
    print(f"{step}: {', '.join(result)}")

def synthesize_detailed():
    state = {
        '[A]': 1, '[B]': 4, '[C]': 7,
        '{A}': 0, '{B}': 0, '{C}': 0,
        '(A)': 0, '(B)': 0, '(C)': 0
    }
    
    print_state(state, "Initial")
    
    while True:
        initial_state = state.copy()
        
        # Try Rule 4 first (combine [C]s)
        while state['[C]'] >= 2:
            state['[C]'] -= 2
            state['{C}'] += 1
            print_state(state, "After Rule 4")
            
        # Try Rule 1
        while state['[A]'] >= 1 and state['[B]'] >= 1 and state['[C]'] >= 1:
            state['[A]'] -= 1
            state['[B]'] -= 1
            state['[C]'] -= 1
            state['{A}'] += 1
            print_state(state, "After Rule 1")
            
        # Continue with other rules...
        # Rule 2
        while state['[A]'] >= 1 and state['[B]'] >= 1:
            state['[A]'] -= 1
            state['[B]'] -= 1
            state['{C}'] += 1
            print_state(state, "After Rule 2")
            
        # Rule 3
        while state['[B]'] >= 1 and state['[C]'] >= 1:
            state['[B]'] -= 1
            state['[C]'] -= 1
            state['{B}'] += 1
            print_state(state, "After Rule 3")
            
        # Rule 5
        while state['{A}'] >= 1 and state['{C}'] >= 1:
            state['{A}'] -= 1
            state['{C}'] -= 1
            state['(A)'] += 1
            state['(B)'] += 1
            print_state(state, "After Rule 5")
            
        # Rule 6
        while state['{B}'] >= 2:
            state['{B}'] -= 2
            state['(C)'] += 1
            print_state(state, "After Rule 6")
            
        if state == initial_state:
            break

synthesize_detailed()