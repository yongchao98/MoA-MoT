import math

def double_factorial(n):
    """Calculates the double factorial of a non-negative integer n."""
    if n < 0:
        raise ValueError("Double factorial is not defined for negative numbers.")
    if n == 0:
        return 1
    result = 1
    # Start from n and decrement by 2 in each step
    for i in range(n, 0, -2):
        result *= i
    return result

# The n-th moment of the conductance g is given by <g^n> = (2n-1)!! / (2n)!!

# 1. Calculate the average conductance <g> (n=1)
n1 = 1
avg_g_num_val = double_factorial(2 * n1 - 1)
avg_g_den_val = double_factorial(2 * n1)
avg_g = avg_g_num_val / avg_g_den_val

print("The problem asks for the ratio <g^4> / <g> for a disordered Majorana wire at the critical point.")
print("The n-th moment of the dimensionless conductance 'g' is given by <g^n> = (2n-1)!! / (2n)!!.")
print("-" * 20)

print(f"First, we calculate the average conductance <g> (for n=1):")
print(f"<g> = ({2*n1 - 1})!! / ({2*n1})!!")
print(f"<g> = {avg_g_num_val} / {avg_g_den_val}")
print(f"<g> = {avg_g}")
print("-" * 20)


# 2. Calculate the fourth moment of the conductance <g^4> (n=4)
n2 = 4
avg_g4_num_val = double_factorial(2 * n2 - 1)
avg_g4_den_val = double_factorial(2 * n2)
avg_g4 = avg_g4_num_val / avg_g4_den_val

print(f"Next, we calculate the fourth moment <g^4> (for n=4):")
print(f"<g^4> = ({2*n2 - 1})!! / ({2*n2})!!")
print(f"<g^4> = {avg_g4_num_val} / {avg_g4_den_val}")
# To show the fraction simplification, we find the greatest common divisor
gcd_val = math.gcd(avg_g4_num_val, avg_g4_den_val)
print(f"<g^4> = {avg_g4_num_val // gcd_val} / {avg_g4_den_val // gcd_val}")
print(f"<g^4> = {avg_g4}")
print("-" * 20)


# 3. Calculate the ratio <g^4> / <g>
ratio = avg_g4 / avg_g
ratio_num = avg_g4_num_val * avg_g_den_val
ratio_den = avg_g4_den_val * avg_g_num_val
gcd_ratio = math.gcd(ratio_num, ratio_den)


print(f"Finally, we calculate the ratio <g^4> / <g>:")
print(f"Ratio = ({avg_g4_num_val} / {avg_g4_den_val}) / ({avg_g_num_val} / {avg_g_den_val})")
print(f"Ratio = {ratio_num} / {ratio_den}")
print(f"Ratio = {ratio_num // gcd_ratio} / {ratio_den // gcd_ratio}")
print(f"The final numerical value of the ratio is: {ratio}")