import sys
import io

# Backup the standard output
stdout_backup = sys.stdout
# Redirect standard output to a string
sys.stdout = io.StringIO()

z = 0.2

# The pdf f_Z(z) is known from derivation to be 6*z^2 - 6*z + 2
a = 6
b = -6
c = 2

term1 = a * (z**2)
term2 = b * z
term3 = c

result = term1 + term2 + term3

print(f"The probability density function (pdf) of Z, f_Z(z), is given by the formula:")
print(f"f_Z(z) = 6z^2 - 6z + 2")
print(f"\nWe need to calculate the value of this function at z = {z}.")
print(f"Substituting z = {z} into the formula:")
print(f"f_Z({z}) = {a}*({z})^2 + ({b})*{z} + {c}")
print(f"\nLet's calculate each term:")
print(f"First term: {a} * ({z})^2 = {a} * {z**2} = {term1:.2f}")
print(f"Second term: {b} * {z} = {term2:.2f}")
print(f"Third term: {c}")
print(f"\nNow, we sum the terms to find the final result:")
print(f"f_Z({z}) = {term1:.2f} + ({term2:.2f}) + {c}")
print(f"f_Z({z}) = {term1 + term2} + {c}")
print(f"f_Z({z}) = {result}")

# Restore the standard output
output = sys.stdout.getvalue()
sys.stdout = stdout_backup
# Print the captured output
print(output)