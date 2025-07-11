import math

# --- User-defined parameters ---
# g is the genus
g = 1
# n is the number of marked points, which is |A|
n = 2
# -----------------------------

print(f"Calculating answers for g={g} and |A|={n}:")
print("---------------------------------------------")

# Part (a): Expression for the minimum number of vertices
# The expression is int(2*g - 2 + |A| > 0)
val_a_interm = 2*g - 2 + n
res_a = int(val_a_interm > 0)
print(f"(a) Minimum number of vertices: ")
print(f"    The calculation is int({2}*{g} - {2} + {n} > 0), which evaluates to int({val_a_interm} > 0), giving the result: {res_a}")
final_a = res_a

# Part (b): Is M_trop_{0,A} always a simplicial fan?
# The theoretical answer is 'yes'.
final_b = "yes"
print(f"\n(b) Is M_trop_{g,A} a simplicial fan for g = 0? \n    The answer is '{final_b}'.")

# Part (c): For g > 0, is it a tropical variety and what is its dimension?
# The theoretical answer to the first question is 'yes'.
final_c1 = "yes"
print(f"\n(c) For g > 0:")
print(f"    Is it a tropical variety? The answer is '{final_c1}'.")
# Its complex dimension is 3g - 3 + |A|, provided the space is non-empty.
if g > 0:
    if val_a_interm > 0:
        res_c2 = 3*g - 3 + n
        print(f"    What is its complex dimension? The calculation is {3}*{g} - {3} + {n}, which gives: {res_c2}")
        final_c2 = res_c2
    else:
        final_c2 = "N/A (space is empty)"
        print(f"    What is its complex dimension? {final_c2}")
else:
    final_c2 = "N/A (question is for g>0)"
    print(f"    What is its complex dimension? {final_c2}")

print("---------------------------------------------")
print("\nFinal answer in the required format for the given values:")
print(f"(a) {final_a}; (b) {final_b}; (c) {final_c1}, {final_c2}")
