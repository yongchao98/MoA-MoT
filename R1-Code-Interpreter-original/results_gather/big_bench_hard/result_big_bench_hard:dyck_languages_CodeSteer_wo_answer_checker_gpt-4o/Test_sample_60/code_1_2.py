# Given input sequence
input_sequence = "{ [ [ [ [ ( ) ] ] ] ]"

# Remove spaces from the input sequence
input_sequence = input_sequence.replace(" ", "")

# Initialize a stack to keep track of opening brackets
stack = []

# Dictionary to match opening and closing brackets
matching_bracket = {')': '(', ']': '[', '}': '{'}

# Iterate through each character in the input sequence
for char in input_sequence:
    if char in '{[(':
        # Push opening brackets onto the stack
        stack.append(char)
        print(f"Pushed {char}, Stack: {stack}")
    elif char in '}])':
        # Pop from stack if the current closing bracket matches the top of the stack
        if stack and stack[-1] == matching_bracket[char]:
            stack.pop()
            print(f"Popped {char}, Stack: {stack}")

# Construct the closing sequence for any unmatched opening brackets
closing_sequence = ''
while stack:
    opening_bracket = stack.pop()
    if opening_bracket == '(':
        closing_sequence += ')'
    elif opening_bracket == '[':
        closing_sequence += ']'
    elif opening_bracket == '{':
        closing_sequence += '}'
    print(f"Closing {opening_bracket}, Stack: {stack}")

# Output the final answer
print(f"<<<{closing_sequence}>>>")