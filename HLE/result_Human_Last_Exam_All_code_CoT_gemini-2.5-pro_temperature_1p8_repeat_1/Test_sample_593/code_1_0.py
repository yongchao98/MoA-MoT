# Symbolic representations for the parameters
t_H = "t_H"
t_G = "t_G"
k = "k"

# The derived tight upper bound for the treewidth of F (t_F)
# t_F is less than or equal to max(t_H, t_G) + k

print("A tight upper bound on the treewidth of F, denoted t_F, is given by the formula:")
print(f"t_F <= max({t_H}, {t_G}) + {k}")