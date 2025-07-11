# n represents the size of the vocabulary, which is 99.
n = 99

# The sum S = sum_{w in V^n} a(w) can be analytically simplified to (n+1)^(n-1).
# Let's calculate the base and exponent of this expression.
base = n + 1
exponent = n - 1

# The result is base^exponent, which is 100^98.
# To express this as a power of 10, we use the property (a^b)^c = a^(b*c).
# Since 100 = 10^2, we have (10^2)^98.
final_base = 10
final_exponent = 2 * exponent

# The final equation is 100^98 = 10^196.
# The following line prints the answer in the requested format (as a power of 10).
# The variables used (base, exponent, final_base, final_exponent) represent
# each number in the final equation.
print(f"{final_base}^{final_exponent}")