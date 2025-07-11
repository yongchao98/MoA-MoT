import collections

# Ranges for the coefficients c_i based on multiplicities in V_n
c1_range = range(3)  # 0, 1, 2
c2_range = range(4)  # 0, 1, 2, 3
c3_range = range(2)  # 0, 1
c4_range = range(2)  # 0, 1

# A set to store the unique coefficient tuples of the independent polynomials
distinct_polynomial_coeffs = set()

# Iterate through all possible combinations of c_1, c_2, c_3, c_4
for c1 in c1_range:
    for c2 in c2_range:
        for c3 in c3_range:
            for c4 in c4_range:
                # The dimension polynomial is p(n) = c1*d1 + c2*d2 + c3*d3 + c4*d4
                # Using the dependency d4 = d1 + d3, we get:
                # p(n) = (c1+c4)*d1 + c2*d2 + (c3+c4)*d3
                # The polynomial is uniquely defined by the coefficients (c1+c4, c2, c3+c4)
                
                c_prime_1 = c1 + c4
                c_prime_2 = c2
                c_prime_3 = c3 + c4
                
                distinct_polynomial_coeffs.add((c_prime_1, c_prime_2, c_prime_3))

# The number of distinct polynomials is the number of unique coefficient tuples.
# This can also be calculated by considering the choices for each new coefficient.
# The choice for c'_2 = c2 is independent of the others.
num_c2_choices = len(c2_range)

# We need to find the number of unique pairs (c'_1, c'_3) = (c1+c4, c3+c4)
c_prime_1_3_pairs = set()
for c1 in c1_range:
    for c3 in c3_range:
        for c4 in c4_range:
            c_prime_1 = c1 + c4
            c_prime_3 = c3 + c4
            c_prime_1_3_pairs.add((c_prime_1, c_prime_3))

num_pair_choices = len(c_prime_1_3_pairs)

# The total number of distinct polynomials is the product of the number of choices.
total_distinct_polynomials = len(distinct_polynomial_coeffs)

print("The number of distinct polynomials is determined by the number of unique coefficient tuples (c1+c4, c2, c3+c4).")
print(f"The number of choices for the second coefficient (from c2) is {num_c2_choices}.")
print(f"The number of unique pairs for the first and third coefficients (from c1, c3, c4) is {num_pair_choices}.")
print(f"The total number of distinct polynomials is the product: {num_c2_choices} * {num_pair_choices} = {total_distinct_polynomials}")
print(f"\nFinal count of distinct polynomials p(n): {total_distinct_polynomials}")