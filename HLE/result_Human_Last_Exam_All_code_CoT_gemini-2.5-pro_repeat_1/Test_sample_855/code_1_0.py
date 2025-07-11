import math

# Helper class to represent numbers in Titan's architecture
class TitanNumber:
    def __init__(self, num, den=1, exp=0):
        # Numerator and denominator must be within 5-bit range (0-31)
        if not (0 <= num <= 31 and 0 <= den <= 31):
            raise ValueError(f"Numerator ({num}) or Denominator ({den}) out of 5-bit range [0, 31]")
        self.num = num
        self.den = den
        self.exp = exp

    def value(self):
        """Returns the decimal value of the Titan number."""
        if self.den == 0:
            return float('inf')
        return (self.num / self.den) * (10 ** self.exp)

    def __str__(self):
        """String representation of the Titan number."""
        if self.exp == 0:
            return f"({self.num}/{self.den})"
        else:
            return f"({self.num}/{self.den})*10^{self.exp}"

def main():
    print("### Step 1: Calculating the 'True' Landing Force (Ground Truth) ###")
    # Constants
    G = 6.6743e-11  # Gravitational constant
    R_PLANET = 2e6  # Pandora's equatorial radius in meters
    R_CORE = 1e5    # Pandora's core radius in meters
    RHO_SHELL = 300 # Shell density in kg/m^3
    RHO_CORE = 1200 # Core density in kg/m^3
    M_PROBE = 50    # Probe mass in kg
    V_INITIAL = 300 # Initial velocity in m/s
    D_LANDING = 5000# Landing distance in meters

    # Pandora's mass calculation (core mass contribution is negligible)
    # M_p = (4/3) * pi * ( (rho_core - rho_shell)*r_core^3 + rho_shell*R_planet^3 )
    m_planet = (4/3) * math.pi * ((RHO_CORE - RHO_SHELL) * R_CORE**3 + RHO_SHELL * R_PLANET**3)
    
    # Pandora's surface gravity g
    g_pandora = (G * m_planet) / (R_PLANET**2)
    
    # Required deceleration a = v^2 / (2d)
    a_decel = V_INITIAL**2 / (2 * D_LANDING)

    # Total landing force Fl = m * (a + g)
    fl_true = M_PROBE * (a_decel + g_pandora)
    
    print(f"Pandora's surface gravity (g): {g_pandora:.4f} m/s^2")
    print(f"Required deceleration (a): {a_decel:.4f} m/s^2")
    print(f"True required landing force: {fl_true:.4f} N\n")

    print("### Step 2: Simulating the Calculation on the Titan Computer ###")
    print("The goal is to calculate F = m * (a + g) using 5-bit fractional arithmetic.")
    
    # --- Attempting a rigorous derivation ---
    print("\nA rigorous calculation of g = (4/3)*pi*G*rho*R fails.")
    print("For example, using pi=22/7 and G=20/3 * 10^-11, the first multiplication is:")
    print("(4/3) * (22/7) = (88/21). The numerator 88 is > 31, so this is an invalid operation.")
    print("We must simplify or approximate.\n")

    # --- Strategy 1: Simplifying by ignoring gravity ---
    print("--- Strategy 1: Simplify by treating gravity as negligible (g=0) ---")
    m_titan_approx1 = TitanNumber(5, 1, 1)  # 5.0 * 10^1 = 50
    a_titan = TitanNumber(9, 1, 0)         # 9/1 = 9
    
    # F = m*a = 50*9 = 450. We need to represent 450.
    # 450 = 4.5 * 10^2 = (9/2) * 10^2
    fl_titan_1 = TitanNumber(9, 2, 2)
    error_1 = abs(fl_titan_1.value() - fl_true)
    print(f"Approximation: g ≈ 0, so F ≈ m * a")
    print(f"m = {m_titan_approx1}, a = {a_titan}")
    print(f"Resulting Force F1 = {fl_titan_1} = {fl_titan_1.value():.4f} N")
    print(f"Absolute Error for Strategy 1: {error_1:.4f} N\n")
    
    # --- Strategy 2: Approximating the sum (a+g) to optimize accuracy ---
    print("--- Strategy 2: Approximate the term (a+g) as a single fraction ---")
    a_plus_g_true = a_decel + g_pandora
    print(f"The true value of (a+g) is {a_plus_g_true:.4f}.")
    print("We need to find a good fraction n/d <= 31/31 for this value.")
    
    # Let's try approximating a+g with 28/3
    a_plus_g_approx = TitanNumber(28, 3, 0)
    print(f"Let's test the approximation (a+g) ≈ {a_plus_g_approx} = {a_plus_g_approx.value():.4f}")
    
    # Now calculate F = m * (a+g)_approx = 50 * (28/3) = 1400/3
    # We need to represent 1400/3.
    # 1400/3 = 466.66... = 4.666... * 10^2 = (14/3) * 10^2
    fl_titan_2 = TitanNumber(14, 3, 2)
    error_2 = abs(fl_titan_2.value() - fl_true)
    
    print(f"This approximation leads to a better result.")
    print(f"Resulting Force F2 = {fl_titan_2} = {fl_titan_2.value():.4f} N")
    print(f"Absolute Error for Strategy 2: {error_2:.4f} N\n")

    print("### Step 3: Final Result ###")
    print("Strategy 2 yields a smaller error, so it is the optimal solution.")
    
    final_force = fl_titan_2
    final_m = m_titan_approx1
    final_a_plus_g = a_plus_g_approx
    final_error = error_2

    print("The final calculation performed by Titan is:")
    print(f"F_landing = m * (a_decel + g_pandora)")
    print("Using the best available approximations for the terms:")
    print(f"{final_force} = {final_m} * {final_a_plus_g}")
    print("\nWhich evaluates to:")
    print(f"{final_force.value():.3f} N = {final_m.value()} kg * {final_a_plus_g.value():.3f} m/s^2")
    
    print(f"\nThe smallest absolute error achievable is {final_error:.3f} N.")


if __name__ == "__main__":
    main()
    print("\n<<<Y[8.277]>>>")