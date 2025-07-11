# The problem reduces to finding the number of regular elements of order 10
# in the reflection group H3. These elements fall into two distinct
# conjugacy classes.

# The group H3 has order 120. Its Coxeter number is 10.

# The first conjugacy class of interest is the Coxeter class.
# Its size is the order of the group divided by the Coxeter number.
num_elements_class1 = 120 // 10

# The second conjugacy class of interest consists of elements conjugate
# to the cube of a Coxeter element. This class is also known to have 12 elements.
num_elements_class2 = 12

# The total number of elements is the sum of the sizes of these two classes.
total_elements = num_elements_class1 + num_elements_class2

# Print the final calculation, showing each number in the equation.
print("The total number of elements is the sum of the sizes of the two regular conjugacy classes of order 10.")
print(f"Number of elements = {num_elements_class1} + {num_elements_class2} = {total_elements}")