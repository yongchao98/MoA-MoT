class TitanOverflowError(Exception):
    """Custom exception for numbers exceeding the 6-bit limit."""
    def __init__(self, value, operation, message="Value exceeds 6-bit limit (0-63) in operation"):
        self.value = value
        self.operation = operation
        self.message = f"{message} '{self.operation}': {self.value}"
        super().__init__(self.message)

class TitanFrac:
    """
    Represents a number in the Titan architecture.
    It's a fraction n/d * 10^e.
    n, d must be within the 6-bit integer range (0-63).
    """
    def __init__(self, n, d=1, e=0, operation="INIT"):
        # Check for overflow upon creation
        for val, name in [(n, "numerator"), (d, "denominator")]:
            if not (0 <= abs(val) <= 63):
                raise TitanOverflowError(val, operation, f"Component '{name}' exceeds 6-bit limit")
        
        self.n = int(n)
        self.d = int(d)
        self.e = int(e)

    def __mul__(self, other):
        # Perform multiplication and check for overflow
        n = self.n * other.n
        d = self.d * other.d
        e = self.e + other.e
        op_str = f"({self}) * ({other})"
        return TitanFrac(n, d, e, operation=op_str)

    def __repr__(self):
        return f"{self.n}/{self.d}e{self.e}"

def solve_pandora_problem():
    """
    Attempts to solve the physics problem using the Titan architecture simulation.
    The function will demonstrate that the calculation is not feasible.
    """
    print("--- Titan Computer Calculation ---")
    print("Goal: Calculate force on a 50kg probe 1km from Pandora black hole's event horizon.")
    print("Pandora: R=2000km, Density=1.2 t/m^3.")
    print("\nStep 1: Define constants and inputs.")
    
    try:
        # Using fractional approximations with components <= 63
        G = TitanFrac(20, 3, -11)
        PI = TitanFrac(22, 7)
        C = TitanFrac(3, 1, 8)
        RHO = TitanFrac(12, 1, 2)
        R = TitanFrac(2, 1, 6)
        FOUR_THIRDS = TitanFrac(4, 3)
        
        print(f"G = {G}, pi = {PI}, c = {C}, rho = {RHO}, R = {R}")
        print("\nStep 2: Attempt to calculate Pandora's Mass M = rho * (4/3) * pi * R^3.")
        
        # This is the first point of failure
        print("Calculating: (rho * 4/3) * pi")
        term1 = RHO * FOUR_THIRDS  # (12e2 * 4/3) = 16e2. This is valid.
        print(f"rho * 4/3 = {term1}")
        
        # The next multiplication causes an overflow
        # (16/1 e2) * (22/7 e0) -> Numerator is 16 * 22 = 352
        # 352 is larger than 63, so this operation fails.
        term2 = term1 * PI
        
    except TitanOverflowError as e:
        print("\n--- CALCULATION FAILED ---")
        print("An intermediate calculation resulted in a number that cannot be represented.")
        print(f"Failure details: {e}")
        print("\nThe calculation of the planet's mass is not possible due to the 6-bit constraint.")
        print("Since the mass (M) cannot be computed, the force F cannot be determined.")
        print("\nFinal Answer: N0")
        print("Each number in the final equation! The equation cannot be formed.")
        print("\n<<<N0>>>")
        return

# Execute the simulation
solve_pandora_problem()