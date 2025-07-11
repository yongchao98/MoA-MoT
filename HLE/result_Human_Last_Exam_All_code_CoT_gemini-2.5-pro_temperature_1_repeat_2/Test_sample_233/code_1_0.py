import math

# The genus of the initial surface Σ is given.
g_Sigma = 10

# Based on the topological argument, the genus of a surface Σ determines the
# maximum possible "linking complexity" of its handles with its boundary.
# A surface of genus g has g handles. In the "worst-case" embedding, each of
# the g handles can be made to link with the boundary component.
# Each such link forces the capping surface S to have a genus of at least 1.
# With 10 handles, the capping surface S might need a genus of up to 10.
g_S_max = 10

# The genus of the final closed surface Σ' is the sum of the genera of the
# two surfaces joined at their boundary.
# g(Σ') = g(Σ) + g(S)
g_prime = g_Sigma + g_S_max

# Output the reasoning and the final calculation, showing each number.
print(f"The problem asks for the smallest integer g such that for any smoothly embedded oriented surface Σ of genus 10 with a single unknotted boundary component, there exists a closed surface Σ' of genus g containing Σ.")
print(f"This can be found by considering the 'worst-case' embedding of Σ.")
print(f"The genus of the final surface, g = g(Σ'), is the sum of the genus of the original surface, g(Σ), and the genus of the 'capping' surface, g(S).")
print(f"The genus of the original surface is g(Σ) = {g_Sigma}.")
print(f"The maximum required genus for the capping surface, for a surface of genus 10, is g(S) = {g_S_max}.")
print(f"Therefore, the final genus is g(Σ') = g(Σ) + g(S).")
print(f"g = {g_Sigma} + {g_S_max} = {g_prime}")
