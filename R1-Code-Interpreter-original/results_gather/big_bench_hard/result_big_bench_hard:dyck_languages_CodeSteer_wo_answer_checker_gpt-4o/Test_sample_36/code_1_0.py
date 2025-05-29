# Initial sequence
sequence = "( ( { }"

# Count the number of opening and closing symbols
open_parentheses = sequence.count('(')
close_parentheses = sequence.count(')')
open_curly_braces = sequence.count('{')
close_curly_braces = sequence.count('}')

# Calculate the number of closing symbols needed
needed_close_parentheses = open_parentheses - close_parentheses
needed_close_curly_braces = open_curly_braces - close_curly_braces

# Complete the sequence
completed_sequence = sequence + ')' * needed_close_parentheses + '}' * needed_close_curly_braces

# Output the final answer
print(f"<<<{completed_sequence}>>>")