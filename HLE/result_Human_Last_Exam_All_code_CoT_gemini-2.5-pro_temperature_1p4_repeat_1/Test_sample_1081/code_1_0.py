# The number of Fq-rational maximal tori of a reductive group of type E8
# over the finite field Fq is equal to the number of conjugacy classes in
# its Weyl group, W(E8). This is because any such group is split,
# and for split groups, this correspondence holds.
# The number of conjugacy classes in W(E8) is a known mathematical constant.

num_conjugacy_classes_W_E8 = 112

# The problem boils down to stating this known number.
# The "equation" is the equality between the number of tori and the number of conjugacy classes.
print(f"Number of Fq-rational maximal tori in G(type E8) = Number of conjugacy classes in W(E8) = {num_conjugacy_classes_W_E8}")
