import math

# Problem parameters based on the interpretation
v = 42  # number of researchers (vertices)
k = 24  # degree of each vertex (number of co-authors)
n_I3 = 2027 # number of tables with no co-authors (independent sets of size 3)

print("Step 1: Setting up the problem.")
print(f"The problem is modeled as a graph with {v} vertices, which is {k}-regular.")
print(f"We are given the number of independent sets of size 3 (tables with no co-authors): n_I3 = {n_I3}.")
print("We need to find the number of triangles (tables where everyone is a co-author), n_K3.\n")

print("Step 2: Using Goodman's theorem for a regular graph.")
print("The formula is: n_K3 + n_I3 = C(v, 3) - v*k*(v-1-k)/2.\n")

print("Step 3: Calculating the components of the formula.")

# Total number of possible 3-person tables
total_triples = math.comb(v, 3)
print(f"Total number of possible tables of 3, C({v}, 3) = {total_triples}.")

# The correction term in the formula
correction_term_numerator = v * k * (v - 1 - k)
correction_term = correction_term_numerator // 2
print(f"The correction term v*k*(v-1-k)/2 = ({v}*{k}*({v}-1-{k}))/2 = {correction_term_numerator}/2 = {correction_term}.")

# The sum of n_K3 and n_I3
total_monochromatic_triangles = total_triples - correction_term
print(f"The sum n_K3 + n_I3 = {total_triples} - {correction_term} = {total_monochromatic_triangles}.\n")

print("Step 4: Solving the final equation for n_K3.")

# Calculate the final answer for n_K3
n_K3 = total_monochromatic_triangles - n_I3

print("The number of tables with all mutual co-authors (n_K3) is calculated as:")
print(f"n_K3 = (n_K3 + n_I3) - n_I3")
print(f"n_K3 = {total_monochromatic_triangles} - {n_I3} = {n_K3}")
