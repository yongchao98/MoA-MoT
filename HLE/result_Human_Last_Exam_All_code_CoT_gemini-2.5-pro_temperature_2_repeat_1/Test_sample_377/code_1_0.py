# The problem reduces to counting the conjugacy classes of the group S = 3^(1+2)_+.
# S is an extraspecial p-group of order p^(1+2n) with p=3 and n=1.

# Order of the group S
order_S = 27

# The center Z(S) of an extraspecial group has order p.
order_Z = 3

# Elements in the center form conjugacy classes of size 1.
num_central_classes = order_Z
print(f"The number of conjugacy classes in the center is: {num_central_classes}")

# Elements outside the center are non-central.
num_non_central_elements = order_S - order_Z
print(f"The number of non-central elements is: {order_S} - {order_Z} = {num_non_central_elements}")

# For any non-central element x, its centralizer C_S(x) is a maximal abelian
# subgroup of order p^(n+1). Here, p=3, n=1, so the order is 3^2=9.
order_centralizer = 9

# The size of a conjugacy class of a non-central element x is |S| / |C_S(x)|.
size_non_central_class = order_S // order_centralizer
print(f"The size of any non-central conjugacy class is: {order_S} / {order_centralizer} = {size_non_central_class}")

# The number of non-central conjugacy classes is the number of non-central
# elements divided by the size of each such class.
num_non_central_classes = num_non_central_elements // size_non_central_class
print(f"The number of non-central conjugacy classes is: {num_non_central_elements} / {size_non_central_class} = {num_non_central_classes}")

# The total number of conjugacy classes is the sum of central and non-central classes.
# This number is equal to the number of blocks of kG.
total_classes = num_central_classes + num_non_central_classes
print(f"\nThe total number of conjugacy classes of S is: {num_central_classes} + {num_non_central_classes} = {total_classes}")
print(f"Therefore, the number of blocks of kG is {total_classes}.")
