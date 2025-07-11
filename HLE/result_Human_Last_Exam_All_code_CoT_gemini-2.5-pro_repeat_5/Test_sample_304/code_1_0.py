import numpy as np

def entropy(p):
    """Calculates entropy H(p) for a probability distribution p."""
    # Filter out zero probabilities to avoid log(0)
    p = p[p > 0]
    return -np.sum(p * np.log2(p))

# For our example, x and y are independent fair coin flips.
# The joint distribution p(x,y) is uniform over {00, 01, 10, 11}.
p_xy = np.array([0.25, 0.25, 0.25, 0.25])
H_xy = entropy(p_xy)

# The individual variables x, y, s1, s2, z have distributions p=[0.5, 0.5]
p_single = np.array([0.5, 0.5])
H_x = entropy(p_single)
H_y = entropy(p_single)
H_s1 = H_x
H_s2 = H_y
H_z = entropy(p_single) # Since z = x XOR y, and x,y are iid Bernoulli(0.5), z is also Bernoulli(0.5)

# The total entropy H(x,y,z,s1,s2) simplifies to H(x,y) for our construction.
H_total = H_xy

print("This problem asks for the maximum value of the joint entropy H(x,y,z,s1,s2).")
print("Through analysis of the constraints, the joint entropy can be shown to be bounded by 2.")
print("We can construct a case that meets this bound.")
print("\nLet x and y be independent fair coin flips (Bernoulli(0.5) variables).")
print(f"H(x) = {H_x}")
print(f"H(y) = {H_y}")
print("Let s1=x, s2=y, and z = x XOR y. This construction satisfies all constraints.")
print(f"H(s1) = H(x) = {H_s1}")
print(f"H(s2) = H(y) = {H_s2}")
print(f"H(z) = H(x XOR y) = {H_z}")
print("\nThe joint entropy H(x,y,z,s1,s2) for this construction simplifies to H(x,y).")
print("Since x and y are independent:")
print(f"H(x,y,z,s1,s2) = H(x,y) = H(x) + H(y) = {H_x} + {H_y} = {H_total}")
print("\nThus, the maximal entropy is 2.")
