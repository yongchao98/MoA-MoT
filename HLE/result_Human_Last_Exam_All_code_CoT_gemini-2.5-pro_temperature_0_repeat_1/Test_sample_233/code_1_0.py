# The genus of the initial surface Sigma.
g_sigma = 10

# The problem asks for the smallest integer g such that a closed surface Sigma'
# of genus g can be formed containing Sigma, regardless of Sigma's embedding.
# A closed surface Sigma' is formed by capping the boundary of Sigma with a surface S.
# The genus of the resulting surface is g_sigma_prime = g_sigma + g_s,
# where g_s is the genus of the capping surface S.

# A key result in topology implies that for any such surface Sigma, it is always
# possible to find a capping surface S that is a disk.
# A disk is a surface of genus 0 with one boundary component.
g_s_min = 0

# Therefore, for any Sigma, we can always construct a closed surface Sigma' of genus:
g_sigma_prime = g_sigma + g_s_min

# The genus of the capping surface g_s must be a non-negative integer, so
# g_sigma_prime can never be less than g_sigma.
# Thus, the smallest such g is the one calculated with the minimal possible g_s.

print(f"The genus of the original surface is {g_sigma}.")
print(f"The minimum possible genus of a valid capping surface is {g_s_min}.")
print("The resulting genus g is the sum of these two values.")
print(f"{g_sigma} + {g_s_min} = {g_sigma_prime}")