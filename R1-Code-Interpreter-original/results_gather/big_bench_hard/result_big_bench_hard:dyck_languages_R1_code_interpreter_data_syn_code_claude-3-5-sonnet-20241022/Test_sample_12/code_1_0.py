def analyze_sequence(input_str):
    # Initialize counters for each bracket type
    count = {'[': 0, ']': 0, '<': 0, '>': 0, '(': 0, ')': 0, '{': 0, '}': 0}
    stack = []
    
    # Process each character
    for char in input_str:
        if char in count:
            count[char] += 1
            if char in '[<({':
                stack.append(char)
            else:
                if stack:
                    # Match closing brackets
                    opening = stack[-1]
                    if ((opening == '[' and char == ']') or
                        (opening == '<' and char == '>') or
                        (opening == '(' and char == ')') or
                        (opening == '{' and char == '}')):
                        stack.pop()
    
    # Generate closing sequence
    result = ''
    while stack:
        bracket = stack.pop()
        if bracket == '[':
            result = ']' + result
        elif bracket == '<':
            result = '>' + result
        elif bracket == '(':
            result = ')' + result
        elif bracket == '{':
            result = '}' + result
    
    print(f"Opening counts: [ = {count['[']}, < = {count['<']}, ( = {count['(']}, {{ = {count['{']}")
    print(f"Closing counts: ] = {count[']']}, > = {count['>']}, ) = {count[')']}, }} = {count['}']}")
    print(f"Required closing sequence: {result}")

input_str = "[ < [ ( ( ) < ( ) > ( { { } } [ [ [ < ( [ ] ) ( ) > ] ] ] { { { { { } } } { { } { < [ [ ] ] > } } { } } } ) ) ] >"
analyze_sequence(input_str)