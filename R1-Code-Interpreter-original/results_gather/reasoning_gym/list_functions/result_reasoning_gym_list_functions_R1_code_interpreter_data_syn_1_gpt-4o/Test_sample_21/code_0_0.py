def find_output(input_list):
    total_sum = sum(input_list)
    for i in range(100):  # Trying different modulo values
        if total_sum % i == 0:
            return i

input_list = [66, 85, 52, 8, 43, 36, 20, 12]
output = find_output(input_list)
print(output)