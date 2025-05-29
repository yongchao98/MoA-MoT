def largest_even_number(numbers):
    even_numbers = [num for num in numbers if num % 2 == 0]
    return max(even_numbers) if even_numbers else None

input_list = [73, 19, 18, 49, 72]
output = largest_even_number(input_list)
print([output])