# The given genus of the initial surface Σ.
genus_Sigma = 10

# The genus of a surface corresponds to its number of handles.
# The problem requires a solution that works for any embedding of Σ.
# We must consider the worst-case scenario where the handles of Σ are
# maximally tangled with respect to its boundary.
# A surface of genus 10 has 10 handles. It is possible to construct
# an embedding where all 10 handles pass through the disk spanned by
# the unknotted boundary circle.
#
# To form a closed surface Σ', we must attach a capping surface, D.
# To avoid intersecting Σ, this capping surface D must have a handle
# for each handle of Σ that it needs to bypass.
# Therefore, the maximum genus of the required capping surface D is 10.
max_genus_capping_surface = 10

# The genus of the final closed surface Σ' is the sum of the genera
# of the two surfaces that are glued together along their boundary.
# g(Σ') = g(Σ) + g(D)
# We use the maximum possible genus for D to find the g that will
# work for all possible Σ.
g = genus_Sigma + max_genus_capping_surface

# Print the final equation and the result.
print(f"g = {genus_Sigma} + {max_genus_capping_surface} = {g}")