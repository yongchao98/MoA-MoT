def find_smallest_odd(numbers):
    odd_numbers = [num for num in numbers if num % 2 != 0]
    if odd_numbers:
        return [min(odd_numbers)]
    return []

input_list = [25, 99, 84, 46, 69, 38, 97, 52, 22]
output = find_smallest_odd(input_list)
print(output)