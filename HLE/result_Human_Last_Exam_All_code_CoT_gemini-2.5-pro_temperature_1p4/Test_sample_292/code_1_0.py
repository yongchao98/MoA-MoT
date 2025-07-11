# The size of the vocabulary.
n = 99

# Based on mathematical derivation, the complex sum simplifies to (n + 1)**(n - 1).
# We calculate the components of this simplified formula.
base = n + 1
exponent = n - 1

# We are asked to write the final answer as a power of 10.
# The base of our result is 100, which is 10**2.
power_of_10_for_base = 2
final_exponent = power_of_10_for_base * exponent

print(f"Given the vocabulary size n = {n}.")
print("The calculation of the sum S = sum_{w} (n+1 - |unique_tokens(w)|)^(-1) simplifies to:")
print(f"S = (n + 1)**(n - 1)")
print(f"Substituting n = {n}, we get:")
print(f"S = ({n} + 1)**({n} - 1)")
print(f"S = {base}**{exponent}")
print("")
print(f"To express this as a power of 10, we first express the base, {base}, as a power of 10:")
print(f"{base} = 10**{power_of_10_for_base}")
print("")
print("Substituting this into the equation for S:")
print(f"S = (10**{power_of_10_for_base})**{exponent}")
print(f"Using the power rule (a**b)**c = a**(b*c), we get:")
print(f"S = 10**({power_of_10_for_base} * {exponent})")
print(f"S = 10**{final_exponent}")