def calculate_expression(y):
    term1 = 114 * 89 * y**2
    term2 = 114 * -57 * y
    term3 = -88 * y * 89 * y**2
    term4 = -88 * y * -57 * y
    
    result = term1 + term2 + term3 + term4
    return result

# Example calculation for y = 1
y_value = 1
result = calculate_expression(y_value)
print(result)