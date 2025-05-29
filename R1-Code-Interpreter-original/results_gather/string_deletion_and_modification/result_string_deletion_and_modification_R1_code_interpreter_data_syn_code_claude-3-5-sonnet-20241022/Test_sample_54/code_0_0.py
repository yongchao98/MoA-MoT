def process_string(s):
    def remove_middle(s):
        if len(s) > 15:
            mid = len(s) // 2
            return s[:mid] + s[mid+1:], True
        return s, False

    def remove_ca_after_first(s):
        if len(s) > 1:
            idx = s[1:].find('ca')
            if idx != -1:
                idx += 1  # Adjust index to account for the slice
                return s[:idx] + s[idx+2:], True
        return s, False

    def remove_bca(s):
        idx = s.find('bca')
        if idx != -1:
            return s[:idx] + s[idx+3:], True
        return s, False

    def replace_abc_with_cab(s):
        idx = s.find('abc')
        if idx != -1:
            return s[:idx] + 'cab' + s[idx+3:], True
        return s, False

    def remove_cc_start(s):
        if s.startswith('cc'):
            return s[2:], True
        return s, False

    def remove_ca_end(s):
        if s.endswith('ca'):
            return s[:-1], True
        return s, False

    operations = [
        remove_middle,
        remove_ca_after_first,
        remove_bca,
        replace_abc_with_cab,
        remove_cc_start,
        remove_ca_end
    ]
    
    current = s
    changes_made = True
    steps = []
    
    while changes_made:
        changes_made = False
        for op in operations:
            new_s, changed = op(current)
            if changed:
                steps.append(f"{current} -> {new_s}")
                current = new_s
                changes_made = True
                break
    
    print("\n".join(steps))
    print("\nFinal string:", current)

# Test the string
process_string("ccbbcabaccacaaccb")