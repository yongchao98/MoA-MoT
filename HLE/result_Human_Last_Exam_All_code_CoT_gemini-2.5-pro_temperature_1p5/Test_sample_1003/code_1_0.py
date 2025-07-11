import math

# The plan is to use the principles of special relativity, specifically the Lorentz
# invariance of the 4-momentum dot product, to relate the angles between stars
# as seen from two different inertial frames. We will set up a system of equations
# based on the given information and solve for the desired ratio.

# Step 1: Establish the key physical formula.
# The transformation of the angle theta_ij between two light sources from a frame 1
# to a frame 2 (primed) is given by the formula:
# 1 - cos(theta'_ij) = (1 - cos(theta_ij)) / (gamma**2 * (1 - v.n_i/c) * (1 - v.n_j/c))
# where v is the relative velocity and n_i are the light direction vectors in frame 1.

# Step 2: Use the information from the first reference frame.
# In frame 1, the angle between any pair of stars is equal, meaning they form a
# regular tetrahedron. The angle 'alpha' between vectors from the center to any
# two vertices of a regular tetrahedron satisfies cos(alpha) = -1/3.
cos_alpha = -1.0 / 3.0
val_1_minus_cos_alpha = 1 - cos_alpha
print("In Frame 1, the stars form a regular tetrahedron.")
print(f"The cosine of the angle between any two stars is {cos_alpha:.4f}.")
print(f"Thus, the term (1 - cos(theta_ij)) is always {val_1_minus_cos_alpha:.4f}.\n")

# Step 3: Set up equations using information from the second reference frame.
# For simplicity, we define a factor D_i = gamma * (1 - v.n_i/c).
# The transformation formula simplifies to: 1 - cos(theta'_ij) = (4/3) / (D_i * D_j).
# We are given the following angles in the second frame:
# - theta'_12 = pi/2, so 1 - cos(theta'_12) = 1.0
# - theta'_13 = 3*pi/4, so 1 - cos(theta'_13) = 1 + 1/sqrt(2)
# - theta'_23 = 3*pi/4, so 1 - cos(theta'_23) = 1 + 1/sqrt(2)
# This gives us a system of equations:
# 1) D1 * D2 = 4/3
# 2) D1 * D3 = (4/3) / (1 + 1/sqrt(2))
# 3) D2 * D3 = (4/3) / (1 + 1/sqrt(2))
print("Using the angles from Frame 2, we establish a system of equations for factors D_i.")
print(f"From S1 and S2: D1 * D2 = {val_1_minus_cos_alpha / 1.0:.4f}")
D1D3 = val_1_minus_cos_alpha / (1 + 1/math.sqrt(2))
print(f"From S1 and S3: D1 * D3 = {D1D3:.4f}")
print("The equation for S2 and S3 is identical to S1 and S3.\n")


# Step 4: Solve for the required ratio.
# The problem asks for the value of (1 - cos(theta'_14)) / (1 - cos(theta'_34)).
# Using our formula, this ratio simplifies:
# Ratio = [ (4/3) / (D1*D4) ] / [ (4/3) / (D3*D4) ] = D3 / D1
# From equations (2) and (3) above, since D1*D3 = D2*D3, we know D1 = D2.
# The ratio D3/D1 can be found by dividing equation (2) by D1^2 (which is D1*D2 = 4/3).
# D3/D1 = (D1*D3)/(D1*D2) = [(4/3)/(1+1/sqrt(2))] / [4/3] = 1/(1+1/sqrt(2))
# Simplifying this expression gives 2 - sqrt(2).
print("The required ratio simplifies to D3 / D1.")
print("Solving the system gives D3 / D1 = 1 / (1 + 1/sqrt(2)), which equals 2 - sqrt(2).\n")

# Step 5: Output the final answer and the numbers in the equation.
final_expression_string = "2 - sqrt(2)"
num_1 = 2.0
num_2_in_sqrt = 2.0
sqrt_of_num_2 = math.sqrt(num_2_in_sqrt)
final_result = num_1 - sqrt_of_num_2

print(f"The final result is calculated from the expression: {final_expression_string}")
print("The numbers in this final equation are:")
print(f"  First number: {num_1}")
print(f"  Number under the square root: {num_2_in_sqrt}")
print("\nFinal calculation:")
print(f"Value = {num_1} - {sqrt_of_num_2}")
print(f"Value = {final_result}")
print("<<<" + f"{final_result}" + ">>>")