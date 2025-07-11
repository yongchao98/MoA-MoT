# The problem reduces to finding the number of conjugacy classes of the group S.
# S is the extraspecial group 3^(1+2)_+, which is of order p^3 with p=3.
p = 3

# The number of conjugacy classes in an extraspecial group of order p^3
# can be calculated by summing the number of central and non-central classes.

# The center has p elements, and each forms its own conjugacy class.
num_central_classes = p

# The number of non-central elements is p^3 - p.
# Each non-central conjugacy class has size p.
# So, the number of non-central classes is (p^3 - p) / p = p^2 - 1.
num_non_central_classes = p**2 - 1

# The total number of conjugacy classes is the sum of these two numbers.
total_classes = num_central_classes + num_non_central_classes

print("The number of blocks of kG is the number of conjugacy classes of S.")
print("This number is the sum of its central and non-central conjugacy classes.")
print(f"Number of central classes = {num_central_classes}")
print(f"Number of non-central classes = {num_non_central_classes}")
print(f"Total number of classes = {num_central_classes} + {num_non_central_classes} = {total_classes}")