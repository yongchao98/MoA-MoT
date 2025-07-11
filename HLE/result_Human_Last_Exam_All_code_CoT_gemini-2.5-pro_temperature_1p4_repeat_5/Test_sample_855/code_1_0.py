import math

# Helper class to simulate Titan's 5-bit fractional arithmetic
class TitanNumber:
    def __init__(self, num, den=1, exp=0):
        if not (0 <= num <= 31 and 0 <= den <= 31):
            raise ValueError(f"Numerator/Denominator {num}/{den} out of 5-bit range (0-31)")
        # Exponent is a signed 5-bit integer, e.g., -16 to 15.
        self.num = num
        self.den = den
        self.exp = exp

    def value(self):
        return (self.num / self.den) * (10 ** self.exp)

    def __str__(self):
        return f"({self.num}/{self.den}) * 10^{self.exp}"

# Define a simplified multiplication that respects Titan's constraints
# In a real scenario, this would involve aggressive simplification/approximation.
# For this simulation, we perform the operation and then create a new TitanNumber.
def titan_multiply(t1, t2):
    # This is a conceptual representation. The actual hardware would need
    # to handle products that exceed the bit limit, likely through approximation.
    new_num = t1.num * t2.num
    new_den = t1.den * t2.den
    new_exp = t1.exp + t2.exp
    
    # Simplify the fraction if possible
    common_divisor = math.gcd(new_num, new_den)
    simplified_num = new_num // common_divisor
    simplified_den = new_den // common_divisor

    # If the simplified result is still out of bounds, an approximation is needed.
    # This is the crucial step where precision might be lost.
    if simplified_num > 31 or simplified_den > 31:
        # Here we would implement an approximation strategy.
        # For this problem, we will show where it fails and manually approximate.
        pass
        
    return simplified_num, simplified_den, new_exp

def main():
    print("Step 1: Define all physical quantities as Titan Numbers.")
    # Probe parameters
    m_probe = TitanNumber(5, 1, 1)  # 50 kg
    v0 = TitanNumber(3, 1, 2)       # 300 m/s
    h = TitanNumber(5, 1, 3)        # 5000 m
    # Pandora parameters
    rho_shell = TitanNumber(3, 1, 2) # 300 kg/m^3
    R = TitanNumber(2, 1, 6)         # 2,000,000 m
    # Constants
    G = TitanNumber(27, 4, -11)      # G ~= 6.75e-11
    pi_approx = TitanNumber(22, 7)   # pi ~= 22/7
    
    print(f"Probe Mass (m): {m_probe} = {m_probe.value()} kg")
    print(f"Initial Velocity (v0): {v0} = {v0.value()} m/s")
    print(f"Altitude (h): {h} = {h.value()} m")
    print(f"Shell Density (rho): {rho_shell} = {rho_shell.value()} kg/m^3")
    print(f"Planet Radius (R): {R} = {R.value()} m")
    print(f"Grav. Constant (G): {G} = {G.value()}")
    print(f"Pi approximation (pi): {pi_approx} = {pi_approx.value()}")
    print("-" * 30)

    print("Step 2: Calculate the kinetic term F_k = (m * v0^2) / (2 * h)")
    v0_squared = TitanNumber(v0.num**2, v0.den**2, v0.exp*2)
    print(f"v0^2 = {v0_squared}")
    
    # Numerator calculation
    # We can simplify m/h = (5e1)/(5e3) = 1e-2
    m_div_h_val = m_probe.value() / h.value()
    # In Titan, we'd simplify (5/1)/(5/1) -> 1/1, and exp 1-3=-2.
    m_div_h = TitanNumber(1, 1, -2)
    print(f"m / h = {m_div_h}")
    
    # Numerator * v0^2
    num_term_val = m_div_h.value() * v0_squared.value()
    num_term_frac_num = m_div_h.num * v0_squared.num
    num_term_frac_den = m_div_h.den * v0_squared.den
    num_term_exp = m_div_h.exp + v0_squared.exp
    print(f"m*v0^2/h = ({num_term_frac_num}/{num_term_frac_den}) * 10^{num_term_exp}")

    # Divide by 2
    Fk_num = num_term_frac_num
    Fk_den = num_term_frac_den * 2
    Fk_exp = num_term_exp
    F_k = TitanNumber(Fk_num, Fk_den, Fk_exp)
    
    print(f"Calculated F_k = {F_k} = {F_k.value()} N")
    print("-" * 30)

    print("Step 3: Attempt to calculate the gravitational term F_g = (4/3)*pi*G*rho*R*m")
    # This involves multiplying six Titan Numbers. Let's trace the fractional parts.
    # Product of fractions: (4/3)*(22/7)*(27/4)*(3/1)*(2/1)*(5/1)
    # The calculation hits a roadblock. E.g., (27/1)*(22/7) -> num=594, which exceeds 31.
    print("Calculation fails: Intermediate products exceed 5-bit limits.")
    print("Example: (4/3)*(3/1) -> (4/1). Then (27/4)*(4/1) -> (27/1).")
    print("Next step: (27/1) * (22/7). Numerator becomes 27*22=594, which is > 31.")
    
    # Approximation is necessary.
    # The actual value of F_g is ~8.4 N.
    # To make the final sum F_k + F_g computable, we approximate F_g.
    print("\nApproximation: The true value of F_g is ~8.4N. We approximate it as 10N.")
    F_g_approx = TitanNumber(10, 1, 0)
    print(f"Approximated F_g = {F_g_approx} = {F_g_approx.value()} N")
    print("-" * 30)
    
    print("Step 4: Calculate the total rocket force F_rocket = F_k + F_g")
    # We need to add F_k=(9/2)*10^2 and F_g=(10/1)*10^0
    Fk_val = F_k.value()
    Fg_val = F_g_approx.value()
    F_rocket_val = Fk_val + Fg_val
    print(f"F_rocket = {Fk_val} N + {Fg_val} N = {F_rocket_val} N")

    # To represent 460 in Titan format, we use scientific notation
    # 460 = 4.6 * 10^2. The fraction for 4.6 is 23/5.
    F_rocket_titan = TitanNumber(23, 5, 2)
    print(f"\nThe final result must be representable. 460 is stored as {F_rocket_titan}.")
    
    print("\nFinal calculation expressed in Titan numbers:")
    print(f"Force = {F_g_approx} + {F_k} = {F_rocket_titan}")
    print(f"Numerically: {F_g_approx.num}/{F_g_approx.den} + ({F_k.num}/{F_k.den})*10^{F_k.exp} = ({F_rocket_titan.num}/{F_rocket_titan.den})*10^{F_rocket_titan.exp}")
    
main()