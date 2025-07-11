from fractions import Fraction

# Step 4: Define the coefficients of the candidate polynomial P(x)
# based on the choice P(1)=1, P(-1)=-1, P(0)=-1.
b_P = Fraction(1 - (-1), 2)
c_P = Fraction(-1)
a_P = Fraction(1 + (-1), 2) - c_P

print(f"The candidate polynomial is P(x) = {a_P}*x^2 + {b_P}*x + {c_P}")

# Step 5: Find the maximum absolute value M of P(x) on [-1, 1].
# The vertex is at x = -b_P / (2*a_P)
vertex_x = -b_P / (2 * a_P)
# Value at the vertex
vertex_val = a_P * vertex_x**2 + b_P * vertex_x + c_P

# The maximum absolute value M is the max of |P(1)|, |P(-1)|, and |P(vertex_x)|
M = max(abs(a_P + b_P + c_P), abs(a_P - b_P + c_P), abs(vertex_val))
print(f"The maximum absolute value of P(x) on [-1, 1] is M = {M}")

# Step 7: Calculate the maximum value for |b| + |c|
# for the scaled polynomial f(x) = P(x)/M.
# The value is (|b_P| + |c_P|) / M.
max_val = (abs(b_P) + abs(c_P)) / M

print("\n--- Final Calculation ---")
print(f"The coefficients of the candidate polynomial are b_P = {b_P} and c_P = {c_P}.")
print(f"The maximum absolute value on [-1, 1] is M = {M}.")
print(f"The maximum value of |b| + |c| is calculated as (|b_P| + |c_P|) / M.")
print(f"|b| + |c| = (|{b_P}| + |{c_P}|) / {M} = {abs(b_P) + abs(c_P)} / {M} = {max_val} = {float(max_val)}")
