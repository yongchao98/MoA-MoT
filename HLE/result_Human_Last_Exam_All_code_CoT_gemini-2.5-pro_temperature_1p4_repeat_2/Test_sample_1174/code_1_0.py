import math

# The derived exponents
n = {
    1: 1,
    2: 1,
    3: -2,
    4: 0,
    5: -2,
    6: -1.5
}

# The expression to calculate is the sum of k * n_k for k from 1 to 6.
total_sum = 0
equation_parts = []

for k in range(1, 7):
    term = k * n[k]
    total_sum += term
    # Prepare parts for printing the full equation
    if n[k] < 0:
        equation_parts.append(f"({k} * ({n[k]}))")
    else:
        equation_parts.append(f"({k} * {n[k]})")

equation_str = " + ".join(equation_parts)
print(f"The values of the exponents are:")
print(f"n1 = {n[1]}")
print(f"n2 = {n[2]}")
print(f"n3 = {n[3]}")
print(f"n4 = {n[4]}")
print(f"n5 = {n[5]}")
print(f"n6 = {n[6]}")
print("\nThe calculation is:")
print(f"{equation_str} = {total_sum}")
