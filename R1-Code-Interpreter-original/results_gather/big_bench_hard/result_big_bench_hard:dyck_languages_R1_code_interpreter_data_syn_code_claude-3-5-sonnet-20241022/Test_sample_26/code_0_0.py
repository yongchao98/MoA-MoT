def is_valid_sequence(seq):
    stack = []
    pairs = {')': '(', '>': '<'}
    
    for char in seq:
        if char in '(<':  # Opening brackets
            stack.append(char)
        elif char in ')>':  # Closing brackets
            if not stack or stack.pop() != pairs[char]:
                return False
    return len(stack) == 0

# Test the sequence with our proposed solution
test_seq = "( ( < > ) )"
print(f"Testing sequence: {test_seq}")
print(f"Is valid: {is_valid_sequence(test_seq)}")

# Let's also verify this is the only valid solution
def generate_possibilities(start):
    base = start.replace(" ", "")
    needed = ""
    stack = []
    # Scan the input to determine what we need
    for char in base:
        if char in '(<':
            stack.append(char)
        elif char in ')>':
            if stack:
                stack.pop()
    # Add needed closing brackets in reverse order
    while stack:
        if stack[-1] == '(':
            needed += ')'
        elif stack[-1] == '<':
            needed += '>'
        stack.pop()
    print(f"The needed closing sequence is: {needed}")
    print(f"Complete sequence should be: {base + needed}")

generate_possibilities("( ( < >")