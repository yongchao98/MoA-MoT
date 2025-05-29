# Correctly evaluate the expression step by step
inner_first = -1 * -6
inner_second = -1 - -8
first_parentheses = inner_first - inner_second

second_inner_first = -7 * 8
second_inner_second = second_inner_first * -4
second_parentheses = 8 + second_inner_second

result = first_parentheses - second_parentheses
print(result)