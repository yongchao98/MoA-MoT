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

def synthesize():
    # Initial state
    state = {
        '[A]': 4, '[B]': 2, '[C]': 2,
        '{A}': 0, '{B}': 0, '{C}': 0,
        '(A)': 0, '(B)': 0, '(C)': 0
    }
    
    while True:
        changes = False
        # Try each method in order
        if apply_method1(state): changes = True
        if apply_method2(state): changes = True
        if apply_method3(state): changes = True
        if apply_method4(state): changes = True
        if apply_method5(state): changes = True
        if apply_method6(state): changes = True
        
        if not changes:
            break
    
    # Create result string
    result = []
    for block, count in state.items():
        if count > 0:
            result.append(f"{count}{block}")
    return " ".join(result)

print(synthesize())