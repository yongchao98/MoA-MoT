def find_missing_brackets(s):
    # Initialize stack and mapping
    stack = []
    opening = "{[(<"
    closing = "}])>"
    brackets_map = {'}': '{', ']': '[', ')': '(', '>': '<'}
    
    # Process the input string
    input_str = s.replace(" ", "")
    
    # First pass: validate and track brackets
    for char in input_str:
        if char in opening:
            stack.append(char)
        elif char in closing:
            if not stack or stack.pop() != brackets_map[char]:
                return "Invalid sequence"
    
    # Generate closing sequence
    missing = ""
    while stack:
        last_open = stack.pop()
        if last_open == '{':
            missing += '}'
        elif last_open == '[':
            missing += ']'
        elif last_open == '(':
            missing += ')'
        elif last_open == '<':
            missing += '>'
    
    print(f"Missing brackets: {missing}")

# Test the input
input_str = "{ < [ < ( < [ ] > ) < > < { } > { < > { } { { < [ < < { } [ { } ] < > > > ] < > > } < [ [ ] ] > } } > ] [ ( ) ] > } [ { ( < ( ( ( [ { } ] ) ) ) > ) } ] [ < ( ) ( ( { ( ) [ { { ( { } ) } } ] } ) ) >"
find_missing_brackets(input_str)