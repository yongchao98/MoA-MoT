def is_prime(n):
    if n <= 1:
        return False
    if n <= 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def digit_sum(n):
    return sum(int(digit) for digit in str(n))

def filter_numbers(numbers):
    result = []
    for number in numbers:
        if is_prime(number) or is_prime(digit_sum(number)):
            result.append(number)
    return result

input_numbers = [85, 56, 1, 71, 26, 85, 21]
output_numbers = filter_numbers(input_numbers)
print(output_numbers)