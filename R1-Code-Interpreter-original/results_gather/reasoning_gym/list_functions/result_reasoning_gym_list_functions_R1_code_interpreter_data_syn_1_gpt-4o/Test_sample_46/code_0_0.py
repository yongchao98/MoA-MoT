def find_median(lst):
    sorted_lst = sorted(lst)
    n = len(sorted_lst)
    if n % 2 == 1:
        return sorted_lst[n // 2]
    else:
        return (sorted_lst[n // 2 - 1] + sorted_lst[n // 2]) / 2

# Test the hypothesis on the given input
input_list = [100, 97, 47, 47, 1]
median = find_median(input_list)
print(median)