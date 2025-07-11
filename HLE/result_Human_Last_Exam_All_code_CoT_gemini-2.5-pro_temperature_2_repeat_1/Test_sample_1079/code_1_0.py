# The reflection group of type H3 is the full symmetry group of the icosahedron.
# Its order can be calculated as the product of its degrees: 2, 6, and 10.
order_H3 = 2 * 6 * 10

# The Coxeter number for H3 is given in the problem statement.
coxeter_number_h = 10

# The number of elements in H3 having a regular eigenvector with an eigenvalue of order 10
# (the Coxeter number) is equal to the size of the conjugacy class of a Coxeter element.
# This is calculated by dividing the order of the group by its Coxeter number.
number_of_elements = order_H3 // coxeter_number_h

# Print the final calculation as requested.
print(f"The number of such elements is calculated by dividing the order of the group H3 by its Coxeter number h.")
print(f"Order of H3 = {order_H3}")
print(f"Coxeter number h = {coxeter_number_h}")
print("The final equation is:")
print(f"{order_H3} / {coxeter_number_h} = {number_of_elements}")