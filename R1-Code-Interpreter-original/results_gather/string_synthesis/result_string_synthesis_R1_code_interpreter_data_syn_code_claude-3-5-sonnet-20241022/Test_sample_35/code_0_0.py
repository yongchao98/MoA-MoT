def try_method_1(state):
    # [A] + [B] + [C] -> {A}
    if state['[A]'] >= 1 and state['[B]'] >= 1 and state['[C]'] >= 1:
        state['[A]'] -= 1
        state['[B]'] -= 1
        state['[C]'] -= 1
        state['{A}'] += 1
        return True
    return False

def try_method_2(state):
    # [A] + [B] -> {C}
    if state['[A]'] >= 1 and state['[B]'] >= 1:
        state['[A]'] -= 1
        state['[B]'] -= 1
        state['{C}'] += 1
        return True
    return False

def try_method_3(state):
    # [B] + [C] -> {B}
    if state['[B]'] >= 1 and state['[C]'] >= 1:
        state['[B]'] -= 1
        state['[C]'] -= 1
        state['{B}'] += 1
        return True
    return False

def try_method_4(state):
    # [C] + [C] -> {C}
    if state['[C]'] >= 2:
        state['[C]'] -= 2
        state['{C}'] += 1
        return True
    return False

def try_method_5(state):
    # {A} + {C} -> (A) + (B)
    if state['{A}'] >= 1 and state['{C}'] >= 1:
        state['{A}'] -= 1
        state['{C}'] -= 1
        state['(A)'] += 1
        state['(B)'] += 1
        return True
    return False

def try_method_6(state):
    # {B} + {B} -> (C)
    if state['{B}'] >= 2:
        state['{B}'] -= 2
        state['(C)'] += 1
        return True
    return False

# Initial state
state = {
    '[A]': 2, '[B]': 4, '[C]': 1,
    '{A}': 0, '{B}': 0, '{C}': 0,
    '(A)': 0, '(B)': 0, '(C)': 0
}

# Keep trying methods until no changes can be made
while True:
    any_change = False
    for method in [try_method_1, try_method_2, try_method_3, try_method_4, try_method_5, try_method_6]:
        if method(state):
            any_change = True
            break
    if not any_change:
        break

# Print final state
result = []
for block_type in ['[A]', '[B]', '[C]', '{A}', '{B}', '{C}', '(A)', '(B)', '(C)']:
    if state[block_type] > 0:
        result.append(f"{state[block_type]} {block_type}")

print(', '.join(result))