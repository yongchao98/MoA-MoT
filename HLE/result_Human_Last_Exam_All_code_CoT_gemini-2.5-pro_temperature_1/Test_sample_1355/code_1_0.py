import math

def double_factorial(n):
    """
    Computes the double factorial n!! recursively.
    For example, 7!! = 7 * 5 * 3 * 1.
    """
    if n < 0:
        raise ValueError("Double factorial is not defined for negative numbers.")
    if n == 0 or n == 1:
        return 1
    return n * double_factorial(n - 2)

# --- Calculation ---

# The k-th moment of the dimensionless conductance 'g' at the critical point
# of a disordered Majorana wire is given by the formula:
# <g^k> = (2k-1)!! / (2k)!!

# 1. Calculate the average conductance <g>, which is the first moment (k=1).
k_avg = 1
avg_g_num = double_factorial(2 * k_avg - 1)
avg_g_den = double_factorial(2 * k_avg)

# 2. Calculate the fourth statistical moment <g^4> (k=4).
k_4th = 4
moment4_g_num_unsimplified = double_factorial(2 * k_4th - 1)
moment4_g_den_unsimplified = double_factorial(2 * k_4th)

# Simplify the fraction for the fourth moment using the greatest common divisor (GCD).
common_divisor_4th = math.gcd(moment4_g_num_unsimplified, moment4_g_den_unsimplified)
moment4_g_num = moment4_g_num_unsimplified // common_divisor_4th
moment4_g_den = moment4_g_den_unsimplified // common_divisor_4th

# 3. Calculate the ratio <g^4> / <g>.
# Ratio = (moment4_g_num / moment4_g_den) / (avg_g_num / avg_g_den)
#       = (moment4_g_num * avg_g_den) / (moment4_g_den * avg_g_num)
ratio_num = moment4_g_num * avg_g_den
ratio_den = moment4_g_den * avg_g_num

# --- Output ---

print("Problem: Find the ratio <g^4> / <g> for a critical Majorana wire.")
print("The k-th moment of conductance 'g' is given by <g^k> = (2k-1)!! / (2k)!!.")
print("-" * 50)

print("Step 1: Calculate the average conductance <g> (for k=1)")
print(f"The equation for the average is: <g> = ({2*k_avg - 1})!! / ({2*k_avg})!!")
print(f"The value of each term is: {avg_g_num} / {avg_g_den}")
print("-" * 50)

print("Step 2: Calculate the fourth moment <g^4> (for k=4)")
print(f"The equation for the fourth moment is: <g^4> = ({2*k_4th - 1})!! / ({2*k_4th})!!")
print(f"The value of each term is: {moment4_g_num_unsimplified} / {moment4_g_den_unsimplified}")
print(f"The simplified fraction is: {moment4_g_num} / {moment4_g_den}")
print("-" * 50)

print("Step 3: Calculate the final ratio <g^4> / <g>")
print(f"The final equation is: ({moment4_g_num}/{moment4_g_den}) / ({avg_g_num}/{avg_g_den}) = {ratio_num}/{ratio_den}")
print(f"The numerical result is {ratio_num/ratio_den}")