import math

def surface_simplicial_volume(g):
    """
    Computes the simplicial volume of an oriented closed surface Sigma_g.
    For g >= 2, ||Sigma_g|| = 4g - 4.
    For g < 2, ||Sigma_g|| = 0.
    """
    if g >= 2:
        return 4 * g - 4
    else:
        return 0

def binomial_coefficient(n, k):
    """
    Computes the binomial coefficient C(n, k).
    """
    if hasattr(math, 'comb'):
        return math.comb(n, k)
    else:
        # For Python versions older than 3.8
        if k < 0 or k > n:
            return 0
        if k == 0 or k == n:
            return 1
        if k > n // 2:
            k = n - k
        
        res = 1
        for i in range(k):
            res = res * (n - i) // (i + 1)
        return res

# Define the genera and dimensions of the surfaces
g1 = 31
g2 = 17
dim1 = 2
dim2 = 2

# Step 1: Compute the simplicial volume of Sigma_31
sv1 = surface_simplicial_volume(g1)

# Step 2: Compute the simplicial volume of Sigma_17
sv2 = surface_simplicial_volume(g2)

# Step 3: Compute the binomial coefficient for the product formula
total_dim = dim1 + dim2
binom_coeff = binomial_coefficient(total_dim, dim1)

# Step 4: Compute the simplicial volume of the product
product_sv = binom_coeff * sv1 * sv2

# Print the results following the formula
print(f"The formula for the simplicial volume of a product of manifolds M x N is:")
print(f"||M x N|| = C(dim(M)+dim(N), dim(M)) * ||M|| * ||N||\n")

print(f"For the surface Sigma_g with genus g >= 2, the simplicial volume is ||Sigma_g|| = 4*g - 4.")
print(f"The simplicial volume of Sigma_{g1} is ||Sigma_{g1}|| = 4*{g1} - 4 = {sv1}")
print(f"The simplicial volume of Sigma_{g2} is ||Sigma_{g2}|| = 4*{g2} - 4 = {sv2}\n")

print(f"The dimensions are dim(Sigma_{g1}) = {dim1} and dim(Sigma_{g2}) = {dim2}.")
print(f"The binomial coefficient is C({dim1}+{dim2}, {dim1}) = C({total_dim}, {dim1}) = {binom_coeff}\n")

print(f"Putting it all together:")
print(f"||Sigma_{g1} x Sigma_{g2}|| = ||Sigma_{g1}|| * ||Sigma_{g2}|| * C({total_dim}, {dim1})")
print(f"= {sv1} * {sv2} * {binom_coeff}")
print(f"= {product_sv}")