# The analysis of the Diophantine equation shows that we need to find the solution
# to x^4 + y^4 + z^4 = w^4 with the smallest maximum value of {x, y, z, w}.
# This is equivalent to finding the solution with the smallest w.

# The smallest known solution in positive integers was discovered by Roger Frye.
# The values for x, y, and z are a permutation of the three bases on the left side.
x = 95800
y = 217519
z = 414560
w = 422481

# The problem asks for the sum of x, y, and z.
sum_of_xyz = x + y + z

# As required, here is the final equation for the sum, with each number outputted.
print(f"{x} + {y} + {z} = {sum_of_xyz}")