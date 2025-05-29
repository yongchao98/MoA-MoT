def find_number_with_highest_second_digit(numbers):
    # Function to extract the second digit
    def second_digit(n):
        return (n // 10) % 10 if n >= 10 else 0
    
    # Find the number with the highest second digit
    max_number = max(numbers, key=lambda x: (second_digit(x), x))
    return [max_number]

# Input list
input_list = [20, 6, 64]
# Find the result
result = find_number_with_highest_second_digit(input_list)
print(result)