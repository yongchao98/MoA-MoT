import math

def compute_gcd(numbers):
    if not all(isinstance(num, int) and num > 0 for num in numbers):
        raise ValueError("All numbers must be positive integers.")
    
    gcd_result = numbers[0]
    for num in numbers[1:]:
        gcd_result = math.gcd(gcd_result, num)
    
    return {"gcd": gcd_result}

# Given input
input_data = {'numbers': [13, 22]}
output = compute_gcd(input_data['numbers'])
print(output)