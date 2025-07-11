import math

# Step 1: Define the entropies of the base variables based on the construction.
# We choose s1 and s2 to be independent Bernoulli(0.5) random variables to maximize their entropy.
# H(X) = -p*log2(p) - (1-p)*log2(1-p). For p=0.5, H(X) = 1.
H_s1 = 1.0
H_s2 = 1.0

# Step 2: According to the analysis, the total joint entropy H(x,y,z,s1,s2)
# simplifies to H(s1, s2) for our construction.
# Since s1 and s2 are independent, H(s1, s2) = H(s1) + H(s2).
H_s1_s2 = H_s1 + H_s2

# Step 3: Print the final calculation, showing each number in the equation.
print(f"The maximal entropy is H(x,y,z,s1,s2) = H(s1,s2).")
print(f"For the optimal construction, s1 and s2 are independent, so H(s1,s2) = H(s1) + H(s2).")
print(f"H(x,y,z,s1,s2) = {H_s1} + {H_s2} = {H_s1_s2}")
