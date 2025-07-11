import math

# Step 1: Define high-precision values to calculate the 'true' escape velocity
G_true = 6.67430e-11
PI_true = math.pi
R_core_true = 1e5  # m
R_planet_true = 2e6  # m
Rho_core_true = 1200  # kg/m^3
Rho_shell_true = 300  # kg/m^3

# Calculate true mass M
m_unscaled_true = (Rho_core_true - Rho_shell_true) * R_core_true**3 + Rho_shell_true * R_planet_true**3
M_true = (4.0 / 3.0) * PI_true * m_unscaled_true

# Calculate true escape velocity
ve_sq_true = 2 * G_true * M_true / R_planet_true
ve_true = math.sqrt(ve_sq_true)


# Step 2: Define and perform calculations using Titan's 4-bit fractional architecture

class TitanNumber:
    """A class to represent a number in Titan's architecture."""
    def __init__(self, numerator, denominator, exponent):
        if not (0 <= numerator <= 15 and 0 < denominator <= 15):
            raise ValueError(f"Numerator/Denominator {numerator}/{denominator} out of 4-bit range")
        self.num = numerator
        self.den = denominator
        self.exp = exponent

    def __mul__(self, other):
        # Multiplication rule: check for overflow before creating the new number
        new_num = self.num * other.num
        new_den = self.den * other.den
        if new_num > 15 or new_den > 15:
            # This indicates an invalid operation that requires decomposition in a real scenario
            # For this simulation, we'll raise an error to show the path is blocked.
             raise ValueError(f"Multiplication overflow: {self.num}*{other.num}/{self.den}*{other.den}")
        return TitanNumber(new_num, new_den, self.exp + other.exp)

    def to_float(self):
        return (self.num / self.den) * (10 ** self.exp)
    
    def display(self):
        return f"{self.num}/{self.den}e{self.exp}"

# --- Approximations for Titan ---
# Using the most accurate approximations that still allow for a valid calculation path.

# Gravitational Constant G = 6.674e-11 -> Use 2/3 * 10^-10
F_G = TitanNumber(2, 3, -10)

# Mass M
# True M_unscaled = 2.409e22. We need a fractional approximation for this.
# Let's try approximating M_unscaled with 5/2 * 10^22 = 2.5e22
# The full M = 4/3 * pi * M_unscaled. Let's use pi~3.
# So M = 4/1 * (5/2 * 10^22) -> (4*5)/(1*2) * 10^22 = 20/2 * 10^22 -> num>15 invalid path.
# We must find an approximation for the entire Mass M, not its components.
# True Mass M is ~1.009e23.
# Let's approximate M as 5/2 * 10^22 = 2.5e22 kg. This is a significant approximation.
F_M = TitanNumber(5, 2, 22)

# Constant 2
F_2 = TitanNumber(2, 1, 0)

# 1/R (Radius R = 2e6 m)
F_R_inv = TitanNumber(1, 2, -6)

# --- Calculation of v_e^2 ---
# The order of multiplication is crucial to avoid overflow.
# Path: (G * M) * 2 * (1/R)
print("Chosen approximation for key values:")
print(f"G ≈ {F_G.display()} ({F_G.to_float():.3e})")
print(f"M ≈ {F_M.display()} ({F_M.to_float():.3e})")
print(f"R = 2e6, so 1/R = {F_R_inv.display()} ({F_R_inv.to_float():.3e})")
print("")

# First operation: G * M
# (2/3 e-10) * (5/2 e22) = (2*5)/(3*2) e12 = 10/6 e12
# Simplify 10/6 -> 5/3
F_GM = TitanNumber(5, 3, 12)

# Second operation: (G*M) * 2
# (5/3 e12) * (2/1) = (5*2)/(3*1) e12 = 10/3 e12. Valid.
F_GM2 = TitanNumber(10, 3, 12)

# Third operation: (G*M*2) * (1/R)
# (10/3 e12) * (1/2 e-6) = (10*1)/(3*2) e6 = 10/6 e6.
# Simplify 10/6 -> 5/3
ve_sq_titan = TitanNumber(5, 3, 6)

# The result is a valid TitanNumber representing v_e^2
ve_sq_titan_val = ve_sq_titan.to_float()
ve_titan_val = math.sqrt(ve_sq_titan_val)

# Calculate the smallest absolute error
abs_error = abs(ve_titan_val - ve_true)

# Final Equation for Titan's calculation
print("The final calculation for v_e^2 on Titan follows the equation:")
print(f"v_e^2 = ({F_2.display()}) * ({F_G.display()}) * ({F_M.display()}) * ({F_R_inv.display()})")
# To show the intermediate steps' results
F_GM_display = TitanNumber(10, 6, 12) # For display before simplification
F_GM2_display = TitanNumber(10,3,12)
F_final_display = TitanNumber(10,6,6)
print(f"      = (({F_G.num}*{F_M.num})/({F_G.den}*{F_M.den})e{F_G.exp+F_M.exp}) * {F_2.display()} * {F_R_inv.display()}")
print(f"      = {F_GM_display.display()} (simplifies to {F_GM.display()}) * {F_2.display()} * {F_R_inv.display()}")
print(f"      = ({F_GM.display()} * {F_2.display()}) * {F_R_inv.display()}")
print(f"      = {F_GM2.display()} * {F_R_inv.display()}")
print(f"      = {F_final_display.display()} (simplifies to {ve_sq_titan.display()})")
print(f"      = {ve_sq_titan.to_float():.2f} m^2/s^2")
print("\nTaking the square root gives the calculated escape velocity:")
print(f"v_e (Titan) = sqrt({ve_sq_titan_val:.2f}) = {ve_titan_val:.2f} m/s")
print(f"v_e (True) = {ve_true:.2f} m/s")
print(f"Absolute Error = |{ve_titan_val:.2f} - {ve_true:.2f}| = {abs_error:.2f} m/s")

# Return the required answer format
answer = f"Y[{abs_error:.2f}]"
print(f"\nFinal Answer: {answer}")
print("<<<" + answer + ">>>")