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

# Initial state
state = {
    '[A]': 3, '[B]': 1, '[C]': 2,
    '{A}': 0, '{B}': 0, '{C}': 0,
    '(A)': 0, '(B)': 0, '(C)': 0
}

methods = [apply_method1, apply_method2, apply_method3, apply_method4, apply_method5, apply_method6]

# Keep applying methods until no changes can be made
while True:
    any_change = False
    for method in methods:
        if method(state):
            any_change = True
    if not any_change:
        break

# Print final state
result = []
for block, count in state.items():
    if count > 0:
        result.append(f"{count} {block}")
print(", ".join(result))