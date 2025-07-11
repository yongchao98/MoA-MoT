from fractions import Fraction

# The order of the group A_5
group_order = 60

# Number of homomorphisms and their conductors for each class
# (unramified, tame C_3, wild C_2, wild V_4)
map_counts = [60, 60, 60, 30]
conductors = [0, 2, 3, 5]

# Calculate the sum of p^(-c) over all homomorphisms
total_sum = 0
for count, c in zip(map_counts, conductors):
    term = count * Fraction(1, 2**c)
    total_sum += term

# Calculate the total mass M = (1/|G|) * sum
mass = Fraction(1, group_order) * total_sum

# Output the equation step-by-step
print("The total mass is calculated as:")
equation_str = f"M(A_5, rho, 2) = (1/{group_order}) * ("
terms_str = []
for i in range(len(map_counts)):
    terms_str.append(f"{map_counts[i]} * 2^(-{conductors[i]})")
equation_str += " + ".join(terms_str) + ")"
print(equation_str)

term_values = []
for i in range(len(map_counts)):
    term_values.append(map_counts[i] * (2**(-conductors[i])))

print(f"= (1/{group_order}) * ({' + '.join(map(str, term_values))})")
print(f"= (1/{group_order}) * {float(total_sum)}")
print(f"= {float(mass)}")
print(f"As a fraction in lowest terms, the total mass is: {mass.numerator}/{mass.denominator}")
