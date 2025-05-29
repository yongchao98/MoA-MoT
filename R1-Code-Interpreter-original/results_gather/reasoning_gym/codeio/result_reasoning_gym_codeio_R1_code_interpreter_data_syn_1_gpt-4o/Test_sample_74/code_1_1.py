def reverse_in_groups(lst, group_size):
    result = []
    for i in range(0, len(lst), group_size):
        result.extend(lst[i:i+group_size][::-1])
    return result

output_list = [8, 61, 68, 72, 3, 47, 64, 25, 9, 34, 37, 84, 95, 58, 56, 36, 60, 96, 29, 89]
group_size = 4

original_list = reverse_in_groups(output_list, group_size)
print(original_list)