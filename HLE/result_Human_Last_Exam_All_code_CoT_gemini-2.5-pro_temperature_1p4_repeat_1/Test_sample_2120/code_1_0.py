# The locations 'z' are the roots of the polynomial originating from the denominator
# of the equation for B(z).
# The polynomial is 4*z**4 - z**3 + z**2 + 1 = 0.
# We identify the coefficients of the polynomial P(z) = a_4*z**4 + a_3*z**3 + a_2*z**2 + a_1*z + a_0.
a_4 = 4
a_3 = -1
a_2 = 1
a_1 = 0
a_0 = 1

# The degree of the polynomial is the number of roots.
n = 4

# According to Vieta's formulas, the sum of the roots of a polynomial is -a_{n-1}/a_n.
# In our case, the sum of roots is -a_3/a_4.
sum_of_roots = -a_3 / a_4

# The average value of the roots is their sum divided by the number of roots.
average_of_roots = sum_of_roots / n

print(f"The polynomial equation for the singularities is: {a_4}z^4 + ({a_3})z^3 + {a_2}z^2 + {a_0} = 0")
print(f"The sum of the roots is calculated as -({a_3})/({a_4}) = {sum_of_roots}")
print(f"The number of roots is {n}")
print(f"The average value of the roots is the sum divided by the number of roots: {sum_of_roots}/{n} = {average_of_roots}")

print("\nThe final calculated average value is:")
print(average_of_roots)