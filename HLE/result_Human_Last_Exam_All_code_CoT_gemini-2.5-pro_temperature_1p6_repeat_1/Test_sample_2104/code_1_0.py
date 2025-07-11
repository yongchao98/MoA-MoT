import math

# Step 1: Define n1 and n2
# Based on the analysis, we exclude the trivial case n=2.
# The 1st and 2nd smallest positive integers satisfying the conditions are 4 and 6.
n1 = 4
n2 = 6

# Step 2: Calculate the parameters of the Hamiltonian's potential
# V(q) = (1/2)q^2 - C * q^k, where H = p^2/2 + V(q)
# H = 1/2 * (p^2 + q^2 - (2/n1)*sqrt((n2-n1)/(n1/2)) * q^(n1/2))
# V(q) = 1/2 * (q^2 - (2/n1)*sqrt((n2-n1)/(n1/2)) * q^(n1/2))
# V(q) = 1/2 * q^2 - (1/n1)*sqrt((n2-n1)/(n1/2)) * q^(n1/2)
# With n1=4, n2=6, the power of q is n1/2=2. The potential is quadratic.
# V(q) = (1/2)kq^2
# Let's calculate the effective spring constant k.
coeff_term_sqrt = math.sqrt((n2 - n1) / (n1 / 2))
coeff = (1 / n1) * coeff_term_sqrt
# The potential V(q) is (1/2 - coeff) * q^2
# So, 1/2 * k = 1/2 - coeff
k = 1 - 2 * coeff

# Step 3: Calculate the period T
# For a simple harmonic oscillator H = p^2/(2m) + (1/2)kq^2
# The mass m = 1
m = 1.0

# Angular frequency omega = sqrt(k/m)
omega = math.sqrt(k / m)

# Period T = 2*pi/omega
T_value = 2 * math.pi / omega

# Step 4: Display the final result
# The problem asks for T((n1-1)/n2)
alpha_num = n1 - 1
alpha_den = n2

# We present the final equation with the numbers plugged in.
# The final result is T(alpha) which is the constant period T_value.
print(f"Let n1 = {n1} and n2 = {n2}.")
print(f"The argument for T is ({n1}-1)/{n2} = {alpha_num}/{alpha_den}.")
print("The effective spring constant k of the resulting Hamiltonian potential is 1/2.")
print("The mass m is 1.")
# The prompt requires showing the numbers in the final equation.
print(f"The final calculation is:")
print(f"T({alpha_num}/{alpha_den}) = 2 * pi / sqrt(k/m) = 2 * pi / sqrt(({k})/({m})) = {T_value}")
