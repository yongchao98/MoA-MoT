# The maximal entropy H(x,y,z,s1,s2) is determined by constructing a
# valid distribution that maximizes the value.

# A known construction uses two independent binary random variables, u and v,
# with H(u) = 1 and H(v) = 1.
# The variables are defined as:
# s1 = u, y = v, x = u XOR v, s2 = u, z = u
# This construction satisfies all the given constraints.

# The total entropy of the system is then the entropy of the base variables:
# H(x,y,z,s1,s2) = H(u,v)
# Since u and v are independent, H(u,v) = H(u) + H(v).
# In our construction H(u) corresponds to H(s1) and H(v) corresponds to H(y).
# All individual entropies are set to their maximum value of 1.

val_s1 = 1.0
val_y = 1.0
max_H = val_s1 + val_y

# The problem asks to output each number in the final equation.
# The equation for the maximal entropy in our construction is H_max = H(s1) + H(y).
print(f"{max_H} = {val_s1} + {val_y}")