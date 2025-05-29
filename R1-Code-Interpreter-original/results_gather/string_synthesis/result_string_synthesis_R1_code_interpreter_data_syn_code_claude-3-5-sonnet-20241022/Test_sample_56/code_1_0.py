def apply_rules(blocks):
    # Deep copy the blocks dictionary for manipulation
    def try_rule1(b):  # [A] + [B] + [C] -> {A}
        if b['[A]'] >= 1 and b['[B]'] >= 1 and b['[C]'] >= 1:
            b['[A]'] -= 1
            b['[B]'] -= 1
            b['[C]'] -= 1
            b['{A}'] += 1
            return True
        return False

    def try_rule2(b):  # [A] + [B] -> {C}
        if b['[A]'] >= 1 and b['[B]'] >= 1:
            b['[A]'] -= 1
            b['[B]'] -= 1
            b['{C}'] += 1
            return True
        return False

    def try_rule3(b):  # [B] + [C] -> {B}
        if b['[B]'] >= 1 and b['[C]'] >= 1:
            b['[B]'] -= 1
            b['[C]'] -= 1
            b['{B}'] += 1
            return True
        return False

    def try_rule4(b):  # [C] + [C] -> {C}
        if b['[C]'] >= 2:
            b['[C]'] -= 2
            b['{C}'] += 1
            return True
        return False

    def try_rule5(b):  # {A} + {C} -> (A) + (B)
        if b['{A}'] >= 1 and b['{C}'] >= 1:
            b['{A}'] -= 1
            b['{C}'] -= 1
            b['(A)'] += 1
            b['(B)'] += 1
            return True
        return False

    def try_rule6(b):  # {B} + {B} -> (C)
        if b['{B}'] >= 2:
            b['{B}'] -= 2
            b['(C)'] += 1
            return True
        return False

    blocks_init = {
        '[A]': 3, '[B]': 4, '[C]': 6,
        '{A}': 0, '{B}': 0, '{C}': 0,
        '(A)': 0, '(B)': 0, '(C)': 0
    }

    # First, try to make as many {A} as possible using rule 1
    while try_rule1(blocks_init):
        pass

    # Then try rule 4 to combine [C]s
    while try_rule4(blocks_init):
        pass

    # Then try rule 3 to make {B}s
    while try_rule3(blocks_init):
        pass

    # Finally try rule 5
    while try_rule5(blocks_init):
        pass

    result = ""
    for block_type in ['[A]', '[B]', '[C]', '{A}', '{B}', '{C}', '(A)', '(B)', '(C)']:
        if blocks_init[block_type] > 0:
            result += f"{blocks_init[block_type]}{block_type} "
    print(result.strip())

apply_rules({})