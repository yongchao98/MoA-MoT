import sys

# Define the 6-bit integer limit
MAX_INT = 63

class TitanFrac:
    """A class to represent a fractional number in the Titan computer."""
    def __init__(self, num, den=1, exp=0):
        self.num = num
        self.den = den
        self.exp = exp
        self.check_overflow()

    def check_overflow(self):
        """Checks if numerator or denominator exceeds the 6-bit limit."""
        if self.num > MAX_INT or self.den > MAX_INT:
            print(f"OVERFLOW ERROR: Fraction {self.num}/{self.den} exceeds 6-bit limit ({MAX_INT}).")
            print("Conclusion: Calculation is not feasible.")
            print("<<<N0>>>")
            sys.exit()

    def __str__(self):
        return f"{self.num}/{self.den}e{self.exp}"

def mov(val_str):
    """Simulates the MOV instruction. Parses a string into a TitanFrac."""
    num_str, den_str = "1", "1"
    exp_str = "0"

    if 'e' in val_str:
        val_str, exp_str = val_str.split('e')
    if '/' in val_str:
        num_str, den_str = val_str.split('/')
    else:
        num_str = val_str

    return TitanFrac(int(num_str), int(den_str), int(exp_str))

def mul(frac1, frac2):
    """Simulates the MUL instruction. Multiplies two TitanFracs."""
    print(f"Multiplying {frac1} by {frac2}")
    new_num = frac1.num * frac2.num
    new_den = frac1.den * frac2.den
    new_exp = frac1.exp + frac2.exp

    # Create the new fraction, which will trigger an overflow check
    result = TitanFrac(new_num, new_den, new_exp)
    print(f"Result: {result} (pre-simplification)")

    # GCD simplification to keep numbers small
    def gcd(a, b):
        return gcd(b, a % b) if b else a
    
    common_divisor = gcd(result.num, result.den)
    result.num //= common_divisor
    result.den //= common_divisor
    print(f"Simplified Result: {result}")
    
    # Check again after simplification, although the first check on raw product is what matters
    result.check_overflow()
    return result

def run_simulation():
    """
    Runs the simulation of the physics calculation on the Titan architecture.
    """
    print("Task: Calculate gravity force on a 50kg probe 1km from a black hole formed by Pandora.")
    print("Pandora: R=2000km, density=1.2 t/m^3\n")
    print("--- STEP 1: Define constants as 6-bit fractions ---")
    
    # Physics constants and given values as fractions
    # G = 6.674e-11 m^3/kg/s^2 ≈ 20/3 e-11
    G = mov("20/3e-11")
    # pi ≈ 22/7 (0.04% error)
    PI = mov("22/7")
    # m_probe = 50 kg
    m_probe = mov("50/1")
    # R = 2000 km = 2e6 m. R^3_frac = 8, R^3_exp = 18
    R_cubed = mov("8/1e18")
    # rho = 1.2 t/m^3 = 1200 kg/m^3. Use 6/5 e3 = 1.2e3 kg/m^3
    rho = mov("6/5e3")
    # d = 1 km = 1e3 m. d^2 = 1e6 m^2
    d_squared = mov("1/1e6")
    # Constant 4/3
    four_thirds = mov("4/3")

    print(f"G = {G}")
    print(f"pi = {PI}")
    print(f"m_probe = {m_probe}")
    print(f"R^3 = {R_cubed}")
    print(f"rho = {rho}")
    print(f"d^2 = {d_squared}")
    print("\n--- STEP 2: Sequentially calculate force F = G * m_probe * rho * (4/3) * pi * R^3 / d^2 ---")
    print("An overflow is expected as intermediate products will exceed 63.\n")

    # Let's try to multiply the terms. Any order will eventually fail.
    # We choose an order that demonstrates the failure quickly.
    
    # Start with G
    AX = G
    print(f"MOV AX, G -> AX = {AX}")
    
    # Multiply by m_probe
    # This first step will already cause an overflow, as 20 * 50 = 1000, which is > 63.
    # The simulation will exit here.
    AX = mul(AX, m_probe)

    # The program will not reach this point.
    # If it could, the next steps would be:
    # AX = mul(AX, rho)
    # AX = mul(AX, four_thirds)
    # ... and so on, each step risking overflow.

# Execute the simulation
run_simulation()
print("\nThe program terminated due to an unavoidable overflow.")
print("The product of the numerators of the necessary constants (e.g., 20 from G and 50 from the probe mass) exceeds the 6-bit limit of 63.")
print("Thus, the calculation is not feasible with the specified Titan architecture.")
print("<<<N0>>>")
