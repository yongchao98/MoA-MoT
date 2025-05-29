def can_apply_method(counts, method_num):
    if method_num == 1:
        return counts['[A]'] >= 1 and counts['[B]'] >= 1 and counts['[C]'] >= 1
    elif method_num == 2:
        return counts['[A]'] >= 1 and counts['[B]'] >= 1
    elif method_num == 3:
        return counts['[B]'] >= 1 and counts['[C]'] >= 1
    elif method_num == 4:
        return counts['[C]'] >= 2
    elif method_num == 5:
        return counts['{A}'] >= 1 and counts['{C}'] >= 1
    elif method_num == 6:
        return counts['{B}'] >= 2
    return False

def apply_method(counts, method_num):
    if method_num == 1:
        counts['[A]'] -= 1
        counts['[B]'] -= 1
        counts['[C]'] -= 1
        counts['{A}'] += 1
    elif method_num == 2:
        counts['[A]'] -= 1
        counts['[B]'] -= 1
        counts['{C}'] += 1
    elif method_num == 3:
        counts['[B]'] -= 1
        counts['[C]'] -= 1
        counts['{B}'] += 1
    elif method_num == 4:
        counts['[C]'] -= 2
        counts['{C}'] += 1
    elif method_num == 5:
        counts['{A}'] -= 1
        counts['{C}'] -= 1
        counts['(A)'] += 1
        counts['(B)'] += 1
    elif method_num == 6:
        counts['{B}'] -= 2
        counts['(C)'] += 1

def synthesize():
    # Initial counts
    counts = {
        '[A]': 3, '[B]': 2, '[C]': 6,
        '{A}': 0, '{B}': 0, '{C}': 0,
        '(A)': 0, '(B)': 0, '(C)': 0
    }
    
    while True:
        changes_made = False
        for method in range(1, 7):
            while can_apply_method(counts, method):
                apply_method(counts, method)
                changes_made = True
        if not changes_made:
            break
    
    # Print final counts
    result = ""
    for block in sorted(counts.keys()):
        if counts[block] > 0:
            result += f"{counts[block]} {block}, "
    print(result[:-2])

synthesize()