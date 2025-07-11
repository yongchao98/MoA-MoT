import math

# Based on the problem description, we have the Johnson scheme J(n,k)
# with n=50, k=5. The graph has vertices as k-subsets of an n-set.
# Adjacency is defined by intersection size of 3.
# The adjacency matrix A is denoted A_i where intersection size is k-i.
# So, intersection size 3 means i = k - 3 = 5 - 3 = 2.
n = 50
k = 5
i = 2

# Step 1: Calculate eigenvalues of the adjacency matrix A (A_2)
# The formula for eigenvalues of A_i in J(n,k) is:
# theta_j = P_i(j) = sum_{u=0 to i} (-1)^u * C(j,u) * C(k-j, i-u) * C(n-k-j, i-u)
# for j = 0, 1, ..., k.
thetas = []
for j in range(k + 1):
    val = 0
    for u in range(i + 1):
        try:
            # Using math.comb for binomial coefficients C(n,k)
            term = ((-1)**u * math.comb(j, u) *
                    math.comb(k - j, i - u) *
                    math.comb(n - k - j, i - u))
            val += term
        except ValueError:
            # math.comb(n, k) raises ValueError if k > n or k < 0, which corresponds to a zero term.
            continue
    thetas.append(val)

# Step 2: Calculate eigenvalues of the Seidel matrix S
# The number of vertices is N = C(n,k).
# The Seidel eigenvalues mu_j are related to the adjacency eigenvalues theta_j.
# For the all-ones eigenvector, mu_0 = N - 1 - 2*theta_0.
# For other eigenvectors, mu_j = -1 - 2*theta_j.
N = math.comb(n, k)
mus = []
mus.append(N - 1 - 2 * thetas[0])
for j in range(1, k + 1):
    mus.append(-1 - 2 * thetas[j])

# Step 3: Compute the LCM of the absolute values of the eigenvalues
abs_mus = [abs(val) for val in mus]

def calculate_lcm_of_list(numbers):
    """Computes the LCM of a list of integers."""
    if not numbers:
        return 0
    result = int(numbers[0])
    for i in range(1, len(numbers)):
        result = (result * int(numbers[i])) // math.gcd(result, int(numbers[i]))
    return result

max_order = calculate_lcm_of_list(abs_mus)

# Step 4: Display the result, showing the components of the calculation.
# This requires the prime power components of the final LCM.

# Prime factorizations of the absolute eigenvalues were pre-analyzed
# as a general factorization algorithm would be too slow for the large numbers involved.
# |mu0| = 2098959 = 3 * 699653
# |mu1| = 11001 = 3 * 19 * 193
# |mu2| = 4905 = 3^2 * 5 * 109
# |mu3| = 1225 = 5^2 * 7^2
# |mu4| = 315 = 3^2 * 5 * 7
# |mu5| = 21 = 3 * 7
factorizations = {
    2098959: {3: 1, 699653: 1},
    11001: {3: 1, 19: 1, 193: 1},
    4905: {3: 2, 5: 1, 109: 1},
    1225: {5: 2, 7: 2},
    315: {3: 2, 5: 1, 7: 1},
    21: {3: 1, 7: 1},
}

lcm_factors = {}
for m in abs_mus:
    current_factors = factorizations.get(m, {})
    for p, e in current_factors.items():
        if p not in lcm_factors or e > lcm_factors[p]:
            lcm_factors[p] = e

# Build the equation string as requested
equation_parts = []
# Sort primes for a consistent, readable output
for p in sorted(lcm_factors.keys()):
    e = lcm_factors[p]
    term = p**e
    equation_parts.append(str(term))

# Final output string construction
final_equation = " * ".join(equation_parts) + f" = {max_order}"

print(final_equation)