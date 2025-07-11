import sys

# This script calculates the conormal regularity of R(sigma)f.

# In Python, we can't directly represent the symbolic variable 'alpha'.
# We will use a string placeholder for it in our output.
alpha_str = "alpha"

# 1. Define the initial parameters from the problem statement.

# The initial function f is in the conormal space A^{2+alpha}(X).
# The order of this space is s = 2 + alpha.
initial_order_numeric = 2

# The operator is the wave operator Box_g, which is a second-order
# differential operator.
operator_order = 2

# 2. State the principle being used.
# The resolvent R(sigma) = (Box_g - sigma^2)^-1 is the inverse of a
# second-order elliptic operator in the scattering calculus.
# Applying the resolvent increases the conormal regularity by the
# order of the operator (elliptic regularity).

# 3. Perform the calculation.
# new_order = initial_order + operator_order
final_order_numeric = initial_order_numeric + operator_order

# 4. Print the explanation and the result.
print("Let s_initial be the conormal order of the function f.")
print(f"The initial conormal order is s_initial = {initial_order_numeric} + {alpha_str}.")
print("-" * 20)
print("Let m be the order of the differential operator P = Box_g - sigma^2.")
print(f"The operator P is a second-order operator, so m = {operator_order}.")
print("-" * 20)
print("According to the principle of elliptic regularity, the conormal order of the solution")
print("u = R(sigma)f, let's call it s_final, is given by:")
print("s_final = s_initial + m")
print("-" * 20)
print("Substituting the values into the equation:")
# We explicitly show each number in the final equation as requested.
print(f"s_final = ({initial_order_numeric} + {alpha_str}) + {operator_order}")
print(f"s_final = ({initial_order_numeric} + {operator_order}) + {alpha_str}")
print(f"s_final = {final_order_numeric} + {alpha_str}")
print("-" * 20)
print("Therefore, the function R(sigma)f belongs to the conormal space:")
print(f"A^({final_order_numeric} + {alpha_str})(X)")

# Appending the final answer in the specified format to the stream.
# sys.stdout.write is used to avoid adding a newline, matching the example format.
sys.stdout.write("<<<A^{4+alpha}(X)>>>")