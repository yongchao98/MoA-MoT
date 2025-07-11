import math

def P(x):
    """
    Calculates the value of the polynomial P(x) for an integer x.
    P(X) = (X^5-1)*(X^5-X)*(X^5-X^2)*(X^5-X^3)*(X^5-X^4)
    """
    if not isinstance(x, int) or x < 0:
        raise ValueError("x must be a non-negative integer")
    if x == 1 or x == 0:
        return 0
    # Use integer arithmetic to avoid precision issues with large numbers
    x5 = x**5
    factors = [
        x5 - 1,
        x5 - x,
        x5 - x**2,
        x5 - x**3,
        x5 - x**4
    ]
    result = 1
    for factor in factors:
        result *= factor
    return result

def gcd_list(numbers):
    """
    Calculates the greatest common divisor of a list of numbers.
    """
    if not numbers:
        return 0
    result = numbers[0]
    for i in range(1, len(numbers)):
        result = math.gcd(result, numbers[i])
    return result

# The limit of the sequence is the fixed divisor of the polynomial P(X).
# We can find this by computing the GCD of the first few non-zero values of P(k).
p_values = []
numbers_to_check = [2, 3, 4, 5]
print("Calculating values of P(k):")
for k in numbers_to_check:
    val = P(k)
    p_values.append(val)
    print(f"P({k}) = {val}")

# Calculate the GCD of these values.
limit_g = gcd_list(p_values)

# Print the final calculation step by step
print("\nCalculating the greatest common divisor:")
current_gcd = p_values[0]
print(f"gcd(P(2)) = {current_gcd}")
for i in range(1, len(p_values)):
    print(f"gcd({current_gcd}, P({numbers_to_check[i]})) = gcd({current_gcd}, {p_values[i]}) = {math.gcd(current_gcd, p_values[i])}")
    current_gcd = math.gcd(current_gcd, p_values[i])
    
print(f"\nThe limit of the sequence (g_n) is {limit_g}")
