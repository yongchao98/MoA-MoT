# The orders of the cyclic groups forming the free product.
n1 = 5
n2 = 8
n3 = 2
k = 3

# The rank 'r' of the kernel K is given by the formula:
# r = 1 - (sum of products of n-1 terms) + (k-1) * (product of all k terms)
sum_of_pairwise_products = n2 * n3 + n1 * n3 + n1 * n2
product_of_all = n1 * n2 * n3
rank = 1 - sum_of_pairwise_products + (k - 1) * product_of_all

print("The fundamental group of Y is the free product G = Z_5 * Z_8 * Z_2.")
print("The first homology group of Y is the direct sum H = Z_5 + Z_8 + Z_2.")
print("The kernel K of the Hurewicz map is the kernel of the abelianization map G -> H.")
print("This kernel K is a free group. Its rank 'r' is calculated as follows:")
print("")
print("Let n1, n2, n3 be the orders of the groups.")
print(f"n1 = {n1}")
print(f"n2 = {n2}")
print(f"n3 = {n3}")
print("")
print("The formula for the rank r is:")
print("r = 1 - (n2*n3 + n1*n3 + n1*n2) + (3-1)*(n1*n2*n3)")
print("")
print("Substituting the values:")
print(f"r = 1 - ({n2}*{n3} + {n1}*{n3} + {n1}*{n2}) + 2*({n1}*{n2}*{n3})")
print(f"r = 1 - ({n2*n3} + {n1*n3} + {n1*n2}) + 2*({product_of_all})")
print(f"r = 1 - ({sum_of_pairwise_products}) + {2 * product_of_all}")
print(f"r = {1 - sum_of_pairwise_products} + {2 * product_of_all}")
print(f"r = {rank}")