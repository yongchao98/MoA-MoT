# Titan architecture constants and constraints
MAX_INT = 31

# Problem parameters represented as fractions
# Height h = 10m
h_n, h_d = 10, 1

# Target distance d = 20m (the only distance for which the calculation is possible)
d_n, d_d = 20, 1

# --- Step 1: Approximate constants for Titan ---

# Gravitational constant g ~ 9.8 m/s^2
# We must use g = 10/1 to ensure intermediate calculations for acceleration
# (a_x = 2*g) do not exceed the 5-bit limit (2*g_n <= 31).
g_n, g_d = 10, 1

# Pi ~ 3.14159
# We must use pi = 3/1 to ensure intermediate calculations for mass
# (m = 3*pi/20) do not exceed the 5-bit limit (3*pi_n <= 31).
pi_n, pi_d = 3, 1

print("Titan Computer Simulation: Calculating the Force")
print("-" * 45)
print(f"Approximations used: pi = {pi_n}/{pi_d}, g = {g_n}/{g_d}")

# --- Step 2: Calculate the mass of the rock (m = 3*pi/20) ---

# m_num = 3 * pi_n
m_interim_n = 3 * pi_n
m_interim_d = 1 * pi_d

# Check for overflow before creating the final fraction
if m_interim_n > MAX_INT:
    print("Error: Mass calculation numerator overflowed.")
else:
    # m = (3*pi) / 20
    m_n = m_interim_n
    m_d = 20
    print(f"1. Mass (m): (3/1 * {pi_n}/{pi_d}) / (20/1) = ({m_n}/{m_d}) kg")

# --- Step 3: Calculate the required horizontal acceleration (a_x = d*g/h) ---

# a_x_num = d_n * g_n
a_x_interim_n = d_n * g_n
a_x_interim_d = d_d * g_d

# Check for overflow. In a real Titan, we'd simplify before multiplying.
# Here we demonstrate simplification: a_x = (20/10) * (10/1) = (20/1)
# which is equivalent to (d/h)*g
a_x_n = (d_n // h_n) * g_n
a_x_d = 1
print(f"2. Acceleration (a_x): ({d_n}/{d_d} * {g_n}/{g_d}) / ({h_n}/{h_d}) = ({a_x_n}/{a_x_d}) m/s^2")

# --- Step 4: Calculate the final force (F = m * a_x) ---

# F = (m_n/m_d) * (a_x_n/a_x_d)
# F_num = m_n * a_x_n = 9 * 20 = 180 (overflow)
# We must simplify before multiplying: (9/20) * (20/1) -> 9/1 * 1/1
final_F_n = m_n * (a_x_n // m_d)
final_F_d = 1

print(f"3. Force (F): ({m_n}/{m_d}) * ({a_x_n}/{a_x_d}) = ({final_F_n}/{final_F_d}) N")
print("-" * 45)
print(f"Final calculated force is {final_F_n}/{final_F_d} N.")
