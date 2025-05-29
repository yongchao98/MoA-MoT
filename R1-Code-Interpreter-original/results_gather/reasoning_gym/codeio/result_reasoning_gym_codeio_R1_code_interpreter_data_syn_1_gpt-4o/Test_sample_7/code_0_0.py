def main_solution(elements, reverse=False):
    if reverse:
        return elements[::-1]
    else:
        return elements

# Test the function with the given output
elements = [82, True, 28.0, True]
reverse = False
output = main_solution(elements, reverse)
print(output)