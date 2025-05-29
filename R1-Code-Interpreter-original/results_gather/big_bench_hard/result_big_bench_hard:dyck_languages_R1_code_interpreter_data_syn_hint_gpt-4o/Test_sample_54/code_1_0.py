def is_balanced(sequence):
    stack = []
    matching_bracket = {')': '(', '}': '{', ']': '[', '>': '<'}
    
    for char in sequence:
        if char in "({[<":
            stack.append(char)
        elif char in ")}]>":
            if stack and stack[-1] == matching_bracket[char]:
                stack.pop()
            else:
                return False
    return not stack

sequence = "(<{{}>{{([({})(([(([({[{}]})]))]))])<[[[[[<>][{[[{((<[]>)[([])]<{[]}>)}}]]}{<(<>)>]]]]]>]]})}}}>"
print(is_balanced(sequence))