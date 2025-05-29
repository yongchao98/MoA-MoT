def is_sorted_ascending(lst):
    return lst == sorted(lst)

input_list = [50, 75, 83]
output = [1] if is_sorted_ascending(input_list) else [0]
print(output)