import math

# The prefactor c_n is a combinatorial term appearing in the virial expansion
# of the equation of state. Based on the analysis of the Mayer cluster expansion,
# this prefactor is given by the formula c_n = -(n-1)/n!.
#
# This script calculates the value of c_n for several small integers n (from 2 to 5)
# and prints the result for each, showing the components of the calculation.

print("The prefactor c_n is part of the virial coefficient B_n, where B_n is decomposed as:")
print("B_n = c_n * Lambda_n + B_n_ex")
print("The derivation from cluster theory shows that the prefactor is c_n = -(n-1) / n!")
print("-" * 30)
print("Here are the values of c_n for small n:")
print("-" * 30)

for n in range(2, 6):
    # Calculate the numerator of the fraction
    numerator = n - 1
    # Calculate the denominator using the factorial function
    denominator = math.factorial(n)
    # Calculate the final value of the prefactor
    c_n = -numerator / denominator

    print(f"For n = {n}:")
    print(f"  The equation is: c_{n} = -({n}-1) / {n}!")
    print(f"  Plugging in the numbers: c_{n} = -{numerator} / {denominator}")
    print(f"  The final value is: c_{n} = {c_n}")
    print("-" * 30)
