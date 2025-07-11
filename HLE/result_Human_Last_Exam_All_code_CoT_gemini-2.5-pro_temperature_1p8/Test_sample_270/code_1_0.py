# The mapping class is phi = (D_a * D_b)^n, with n=9.
n = 9

# We use the known relation from the mapping class group of a punctured torus:
# (D_a * D_b)^p = D_delta^q
# where p=6 and q=-1.
p = 6
q = -1

# We want to find the Fractional Dehn Twist Coefficient (FDTC) of phi.
# The definition of FDTC is c = m/k, where phi^k = D_delta^m.

# We can raise phi to the power of p:
# phi^p = ((D_a * D_b)^n)^p = ((D_a * D_b)^p)^n = (D_delta^q)^n = D_delta^(q*n)
# So, we have found our k and m:
k = p
m = q * n

# Now we calculate the coefficient c = m/k.
c = m / k

print("The problem is to find the fractional Dehn twist coefficient 'c' of the mapping class phi = (D_a * D_b)^9.")
print("This coefficient is defined by the relation phi^k = D_delta^m, where c = m/k.")
print("Using the known algebraic relation (D_a * D_b)^6 = D_delta^(-1), we find a suitable power k.")
print("By setting k=6, we get phi^6 = ( (D_a * D_b)^9 )^6 = ( (D_a * D_b)^6 )^9 = (D_delta^(-1))^9 = D_delta^(-9).")
print(f"So, k={k} and m={m}.")
print("The final equation for the coefficient is:")
print(f"c = {m} / {k}")
print(f"The result is: {c}")

# Alternatively, from the derivation in the text:
# We found that phi^2 = D_delta^(-3).
k_alt = 2
m_alt = -3
c_alt = m_alt / k_alt
# The ratio m/k is the same: -9/6 = -3/2.
