import numpy as np

# Step 1: Define the vectors and problem parameters.
# The three types of vectors
v_a = np.array([1, 0])
v_b = np.array([0.5, np.sqrt(3)/2])
v_c = np.array([-0.5, np.sqrt(3)/2])

print("--- Step 1: Problem Setup ---")
print("The problem involves a sum S of n = 6k random vectors.")
print(f"The three types of vectors are:")
print(f"  - v_a = {v_a}")
print(f"  - v_b = {v_b}")
print(f"  - v_c = {v_c}")
print("Each vector is chosen with a random sign (+1 or -1).\n")

# Step 2: Calculate the covariance matrix for the CLT.
# The random vector sum is S = sum(epsilon_i * v_i).
# By the Central Limit Theorem, the distribution of S/sqrt(n) converges to a 
# Normal distribution N(0, Sigma_prime).
# The covariance matrix is Sigma_prime = (1/n) * sum(v_i * v_i^T for i=1 to n).

# Calculate v * v^T for each vector type (as an outer product)
cov_a = np.outer(v_a, v_a)
cov_b = np.outer(v_b, v_b)
cov_c = np.outer(v_c, v_c)

# The sum of these individual covariance matrices for the three unique vectors
sum_cov_indiv = cov_a + cov_b + cov_c

# The total sum of covariance matrices is Sum_Cov = 2k * (cov_a + cov_b + cov_c).
# Since n = 6k, we have 2k = n/3.
# So, Sigma_prime = (1/n) * (n/3) * sum_cov_indiv = (1/3) * sum_cov_indiv
Sigma_prime = (1/3) * sum_cov_indiv

print("--- Step 2: Covariance Matrix for the Limiting Distribution ---")
print("The sum of outer products (v * v^T) for the three unique vectors is:")
print(f"{sum_cov_indiv}")
print("This is exactly the matrix [[1.5, 0], [0, 1.5]] = (3/2) * I.\n")
print(f"The covariance matrix for S/sqrt(n) is Sigma' = (1/3) * (3/2) * I:")
print(f"{Sigma_prime}")
print("This is exactly the matrix [[0.5, 0], [0, 0.5]] = (1/2) * I.\n")

# Step 3: Probability Density Function (PDF)
# The vector Z = S/sqrt(n) follows N(0, Sigma_prime).
# The PDF is f(z) = (1 / (2*pi*sqrt(det(Sigma')))) * exp(-0.5 * z^T * inv(Sigma') * z)
det_Sigma_prime = np.linalg.det(Sigma_prime)
# Analytically, inv(Sigma') is [[2, 0], [0, 2]].

print("--- Step 3: Probability Density of the Limiting Distribution ---")
print("The distribution of Z = S/sqrt(n) approaches a Gaussian N(0, Sigma').")
print(f"The PDF is f(z1, z2) = (1 / (2*pi*sqrt({det_Sigma_prime:.1f}))) * exp(- (z1^2 + z2^2))")
print("Which simplifies to: f(z1, z2) = (1 / pi) * exp(-(z1^2 + z2^2)).\n")


# Step 4: Connecting P(n) to the Limiting Distribution
# We want to find P(n) = P(||S||^2 <= 2).
# This is equivalent to P(||sqrt(n) * Z||^2 <= 2), which simplifies to P(||Z||^2 <= 2/n).
# The probability is found by integrating the PDF f(z1, z2) over a disk of radius sqrt(2/n).
# The exact integral gives P(n) = 1 - exp(-2/n).

print("--- Step 4: The Probability P(n) ---")
print("The probability P(n) = P(||S||^2 <= 2) becomes P(||Z||^2 <= 2/n).")
print("By integrating the Gaussian PDF over this small disk, we find the exact expression for large n:")
print("P(n) = 1 - exp(-2/n)\n")


# Step 5: The Limit
# We need to find the limit of n * P(n) as n -> infinity.
# The expression is: lim_{n->inf} n * (1 - exp(-2/n))
# To see how this limit is evaluated, we can use the Taylor expansion for exp(x)
# around x=0, which is e^x ~ 1 + x.
# Let x = -2/n. For large n, x is close to 0.
# So, exp(-2/n) ~ 1 - 2/n.
# Then, 1 - exp(-2/n) ~ 1 - (1 - 2/n) = 2/n.
# Therefore, n * P(n) ~ n * (2/n) = 2.

print("--- Step 5: The Final Limit ---")
print("We want to compute the limit of n * P(n) as n goes to infinity.")
print("The expression for the limit is: L = lim_{n->inf} [n * (1 - exp(-2/n))]")
print("Let's break down the expression to highlight its components:")

var_n = 'n'
one = 1
minus_two = -2
final_limit = 2.0

print(f"  - The term P({var_n}) is ( {one} - exp({minus_two}/{var_n}) ).")
print(f"  - We are looking for the limit of the product: {var_n} * P({var_n}).")
print(f"As explained, for large n, P({var_n}) approximates to 2/{var_n}.")
print(f"So, the product {var_n} * P({var_n}) approaches {var_n} * (2/{var_n}) = 2.")
print(f"\nThe exact value of the limit is {final_limit}.")
