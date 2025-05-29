def find_first_repeated_index(lst):
    seen = {}
    for index, value in enumerate(lst):
        if value in seen:
            return seen[value]
        seen[value] = index
    return 0

input_list = [5, 38, 5, 5, 68, 71, 71, 3, 100, 5, 21, 57]
output = find_first_repeated_index(input_list)
print([output])