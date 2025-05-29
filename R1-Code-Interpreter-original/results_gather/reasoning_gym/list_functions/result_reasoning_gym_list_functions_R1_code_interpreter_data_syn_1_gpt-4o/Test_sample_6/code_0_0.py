def largest_multiple_of_11(numbers):
    multiples_of_11 = [num for num in numbers if num % 11 == 0]
    return max(multiples_of_11) if multiples_of_11 else None

input_list = [49, 48, 93]
output = largest_multiple_of_11(input_list)
print(output)