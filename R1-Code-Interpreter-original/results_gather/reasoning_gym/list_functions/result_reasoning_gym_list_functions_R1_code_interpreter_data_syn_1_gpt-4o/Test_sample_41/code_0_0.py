def calculate_output(input_list):
    even_count = sum(1 for x in input_list if x % 2 == 0)
    odd_count = len(input_list) - even_count
    if even_count == 0 or odd_count == 0:
        return 0
    else:
        return even_count

input_list = [64, 71, 56, 24, 8, 34, 73, 77, 26]
output = calculate_output(input_list)
print(output)