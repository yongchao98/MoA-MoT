import math

def count_isotropic_vectors(n, d_eff):
    """
    Counts the number of non-zero binary vectors of length n whose
    Hamming weight k is a multiple of d_eff. This corresponds to the number
    of d-primitive, d-isotropic vectors mod d, where d=2 and d_eff=ed=4.
    """
    total = 0
    terms = []
    values = []
    
    # The weight k must be non-zero (for primitivity) and a multiple of d_eff.
    # The loop starts at d_eff because k=0 (the zero vector) is not primitive.
    for k in range(d_eff, n + 1, d_eff):
        term = math.comb(n, k)
        terms.append(f"C({n},{k})")
        values.append(str(term))
        total += term

    print(f"The calculation for n={n} requires summing C(n,k) for k > 0 and k divisible by {d_eff}:")
    equation = " + ".join(terms)
    value_equation = " + ".join(values)
    print(f"{equation} = {value_equation} = {total}")
    return total

# --- Solve Part (a) ---
print("(a) To find the number of distinct 2-isotropic vectors for n=16, we count non-zero vectors in (Z/2Z)^16 whose weight is a multiple of 4.")
val_a = count_isotropic_vectors(n=16, d_eff=4)
print("-" * 20)

# --- Solve Part (b) ---
print("(b) To determine if a 3-isotropic vector x in Z^8 can result in an even 3-neighbor, we analyze the conditions.")
print("The conditions are contradictory. A 3-isotropic vector in Z^8 must have sum of squares divisible by 3. For the neighbor to be even, all components of x must be non-divisible by 3, which implies the sum of squares is 8 = 2 (mod 3). This is a contradiction.")
val_b = "no"
print(f"Therefore, the answer is: {val_b}")
print("-" * 20)

# --- Solve Part (c) ---
print("(c) To find the number of unimodular 2-neighbors of Z^12, we count the number of 2-primitive, 2-isotropic vectors modulo 2.")
print("This is equivalent to counting non-zero vectors in (Z/2Z)^12 whose weight is a multiple of 4.")
val_c = count_isotropic_vectors(n=12, d_eff=4)
print("-" * 20)

# --- Final Answer ---
# The final answer is formatted as requested by the user.
final_answer = f"(a) {val_a}; (b) {val_b}; (c) {val_c}"
print(f"\nFinal Answer: {final_answer}")