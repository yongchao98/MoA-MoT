def apply_method1(state):
    # [A] + [B] + [C] -> {A}
    if state['[A]'] >= 1 and state['[B]'] >= 1 and state['[C]'] >= 1:
        state['[A]'] -= 1
        state['[B]'] -= 1
        state['[C]'] -= 1
        state['{A}'] += 1
        return True
    return False

def apply_method2(state):
    # [A] + [B] -> {C}
    if state['[A]'] >= 1 and state['[B]'] >= 1:
        state['[A]'] -= 1
        state['[B]'] -= 1
        state['{C}'] += 1
        return True
    return False

def apply_method3(state):
    # [B] + [C] -> {B}
    if state['[B]'] >= 1 and state['[C]'] >= 1:
        state['[B]'] -= 1
        state['[C]'] -= 1
        state['{B}'] += 1
        return True
    return False

def apply_method4(state):
    # [C] + [C] -> {C}
    if state['[C]'] >= 2:
        state['[C]'] -= 2
        state['{C}'] += 1
        return True
    return False

def apply_method5(state):
    # {A} + {C} -> (A) + (B)
    if state['{A}'] >= 1 and state['{C}'] >= 1:
        state['{A}'] -= 1
        state['{C}'] -= 1
        state['(A)'] += 1
        state['(B)'] += 1
        return True
    return False

def apply_method6(state):
    # {B} + {B} -> (C)
    if state['{B}'] >= 2:
        state['{B}'] -= 2
        state['(C)'] += 1
        return True
    return False

def simulate():
    # Initial state
    state = {
        '[A]': 4, '[B]': 2, '[C]': 2,
        '{A}': 0, '{B}': 0, '{C}': 0,
        '(A)': 0, '(B)': 0, '(C)': 0
    }
    
    methods = [apply_method1, apply_method2, apply_method3, 
              apply_method4, apply_method5, apply_method6]
    
    while True:
        changes_made = False
        for method in methods:
            if method(state):
                changes_made = True
                break
        if not changes_made:
            break
    
    # Print final state
    result = []
    for block, count in sorted(state.items()):
        if count > 0:
            result.append(f"{block}: {count}")
    print(", ".join(result))

simulate()