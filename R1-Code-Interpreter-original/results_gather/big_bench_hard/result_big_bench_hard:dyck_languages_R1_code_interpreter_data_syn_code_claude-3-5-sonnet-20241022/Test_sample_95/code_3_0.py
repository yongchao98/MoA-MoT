def validate_final_sequence(s):
    pairs = {
        '(': ')',
        '[': ']',
        '{': '}',
        '<': '>'
    }
    stack = []
    sequence = s.replace(" ", "")
    
    # Count brackets
    counts = {char: 0 for char in '()[]{}<>'}
    for char in sequence:
        counts[char] += 1
    
    # Check matching counts
    for open_br, close_br in pairs.items():
        if counts[open_br] != counts[close_br]:
            print(f"Mismatch in counts: {open_br}={counts[open_br]}, {close_br}={counts[close_br]}")
            return False
    
    # Validate sequence
    for i, char in enumerate(sequence):
        if char in '([{<':
            stack.append(char)
        elif char in ')]}>':
            if not stack:
                print(f"Error at position {i}: No matching opening bracket for {char}")
                return False
            if pairs[stack[-1]] != char:
                print(f"Error at position {i}: Expected {pairs[stack[-1]]}, found {char}")
                return False
            stack.pop()
    
    if stack:
        print(f"Unclosed brackets remaining: {''.join(stack)}")
        return False
    
    print("Sequence is valid!")
    print(f"Total pairs: {sum(counts[char] for char in '([{<')} of each type")
    return True

# Test the sequence
sequence = "((([<[<[([{<(<(()){}[]>)>}]){<>}][<[{(<<{{[<{[[<>[]]]}>]{{}}}}>>)}]>{<>{}}]>([])]>])()))"
validate_final_sequence(sequence)