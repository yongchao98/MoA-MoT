import math

def F_func(u, u_bar):
  """Calculates F(u, u_bar) = u(1-u)^2 * e^(-u_bar)"""
  return u * (1 - u)**2 * math.exp(-u_bar)

def F11_func(u, u_bar):
  """Calculates F_11(u, u_bar) = (-4+6u) * e^(-u_bar)"""
  return (-4 + 6 * u) * math.exp(-u_bar)

def A_func(u, u_bar):
  """Calculates A(u, u_bar) = (F_111*F_2 + F_1*F_11)"""
  p_u = 12 * u**3 - 24 * u**2 + 16 * u - 4
  return p_u * math.exp(-2 * u_bar)

# Point of maximum found by analytical investigation of the expression for E
u = 0.0
u_bar = 0.0
u_star = 1.0
u_star_bar = 0.0

# Calculate the components of the expression E at the maximum
F_val = F_func(u, u_bar)
F_star_val = F_func(u_star, u_star_bar)

delta_u = u_star - u
delta_F = F_star_val - F_val

F11_val = F11_func(u, u_bar)
A_val = A_func(u, u_bar)

# The expression to maximize is E = F_11 * delta_F - A * delta_u
max_E = F11_val * delta_F - A_val * delta_u

# Output the equation with the numbers plugged in
print(f"The expression to maximize is E = F_11 * (F(u_*, u_bar_*) - F(u, u_bar)) - A * (u_* - u)")
print("The maximum is found at u=0, u_bar=0, u_*=1, u_bar_*=0.")
print("\nCalculating each term at this point:")
print(f"F_11(0, 0) = {F11_val:.4f}")
print(f"F(1, 0) = {F_star_val:.4f}")
print(f"F(0, 0) = {F_val:.4f}")
print(f"Delta_F = F(1, 0) - F(0, 0) = {delta_F:.4f}")
print(f"A(0, 0) = {A_val:.4f}")
print(f"u_* = {u_star}")
print(f"u = {u}")
print(f"Delta_u = u_* - u = {delta_u:.4f}")
print("\nFinal calculation:")
print(f"E_max = ({F11_val:.4f}) * ({delta_F:.4f}) - ({A_val:.4f}) * ({delta_u:.4f})")
print(f"E_max = {F11_val * delta_F} - {A_val * delta_u}")
print(f"E_max = {max_E:.4f}")

<<<4>>>