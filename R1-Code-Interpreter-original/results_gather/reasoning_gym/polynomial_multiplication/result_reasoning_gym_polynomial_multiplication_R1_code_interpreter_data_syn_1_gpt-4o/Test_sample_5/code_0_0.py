def compute_expression(z):
    term1 = 6930 * z**4
    term2 = 11663 * z**3
    term3 = 4902 * z**2
    term4 = -693 * z
    term5 = -602
    result = term1 + term2 + term3 + term4 + term5
    return result

# Example: Compute the expression for z = 1
z_value = 1
result = compute_expression(z_value)
print(result)