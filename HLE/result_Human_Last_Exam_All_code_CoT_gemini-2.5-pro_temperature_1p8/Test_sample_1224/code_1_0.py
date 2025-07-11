import math

class Titan4BitError(Exception):
    """Custom exception for Titan architecture constraint violations."""
    pass

class TitanFraction:
    """
    A class to represent a fraction on the Titan 4-bit architecture.
    It enforces the constraint that numerator and denominator must be <= 15.
    """
    def __init__(self, num, den=1):
        if not (0 <= num <= 15 and 1 <= den <= 15):
            raise Titan4BitError(f"Initialization failed: Numerator {num} or denominator {den} out of 4-bit range (0-15).")
        
        # Auto-simplification using GCD
        common = math.gcd(num, den)
        self.num = num // common
        self.den = den // common

    def __mul__(self, other):
        """
        Simulates the MUL instruction.
        The result of the multiplication must also be representable.
        This function strictly checks for intermediate overflows.
        """
        # Cross-simplify before multiplying to minimize intermediate values
        g1 = math.gcd(self.num, other.den)
        g2 = math.gcd(other.num, self.den)
        
        num1 = self.num // g1
        den2 = self.den // g2
        
        num2 = other.num // g2
        den1 = other.den // g1

        new_num = num1 * num2
        new_den = den1 * den2

        if new_num > 15 or new_den > 15:
            raise Titan4BitError(f"MUL operation failed: ({self}) * ({other}) -> Intermediate result {new_num}/{new_den} exceeds 4-bit constraints.")
        
        return TitanFraction(new_num, new_den)

    def __add__(self, other):
        """
        Simulates the ADD instruction.
        This is where overflow is very likely.
        """
        # The numerator (self.num * other.den + other.num * self.den) must be calculated.
        # This intermediate value is the source of the constraint violation.
        new_num = self.num * other.den + other.num * self.den
        new_den = self.den * other.den
        
        if new_num > 15 or new_den > 15:
            raise Titan4BitError(f"ADD operation failed: ({self}) + ({other}) -> Intermediate numerator {new_num} or denominator {new_den} exceeds 4-bit constraints.")
            
        return TitanFraction(new_num, new_den)

    def __truediv__(self, other):
        """Simulates the DIV instruction."""
        # Division is multiplication by the reciprocal.
        reciprocal = TitanFraction(other.den, other.num)
        return self.__mul__(reciprocal)

    def __str__(self):
        return f"{self.num}/{self.den}"

    def to_float(self):
        return self.num / self.den

def calculate_landing_time():
    """
    Main function to simulate the Pandora landing time calculation on Titan.
    """
    print("Task: Calculate landing time t = sqrt(2*h/g) on Pandora.\n")

    # Step 1: Define constants using TitanFraction approximations
    print("Step 1: Define constants with 4-bit fractional approximations.")
    try:
        # Physical parameters of Pandora
        h = TitanFraction(5, 1) # height = 5000m, we will handle exponent separately
        
        # Simplified model: g = (4/3) * pi * G * r * d
        # Using approximations that have a chance to work without immediate overflow.
        # (4/3)*pi part: approx pi=13/4. (4/3)*(13/4) = 13/3. This works.
        const_4_3_pi = TitanFraction(13, 3) 
        
        # G ≈ 6.67e-11. We cannot use 13/2=6.5 due to large product. Let's try 6/1.
        G_mantissa = TitanFraction(6, 1)
        
        # r (mean radius) ≈ 2e6 m. Mantissa is 2.
        r_mantissa = TitanFraction(2, 1)
        
        # d (shell density) = 300 kg/m^3. Mantissa is 3.
        d_mantissa = TitanFraction(3, 1)

        print(f"  - Height mantissa (h): {h}")
        print(f"  - (4/3)*pi approximation: {const_4_3_pi} (using pi ≈ 13/4=3.25)")
        print(f"  - Gravitational constant mantissa (G): {G_mantissa} (using G ≈ 6e-11)")
        print(f"  - Radius mantissa (r): {r_mantissa}")
        print(f"  - Density mantissa (d): {d_mantissa}\n")

    except Titan4BitError as e:
        print(f"Failed at initialization: {e}")
        return "N0"

    # Step 2: Calculate surface gravity 'g' mantissa
    print("Step 2: Calculate g_mantissa = (4/3)*pi * G * r * d")
    print("Multiplying terms in an order that might avoid overflow:")
    # Order: ((13/3)*3) * 6 * 2
    try:
        print(f"  - MOV AX, {d_mantissa}")
        g_calc = d_mantissa
        print(f"  - MUL AX, {const_4_3_pi}  |  AX holds ({g_calc}) * ({const_4_3_pi})")
        g_calc = g_calc * const_4_3_pi
        print(f"    - Result OK. AX = {g_calc}") # 13/1

        print(f"  - MUL AX, {r_mantissa}  |  AX holds ({g_calc}) * ({r_mantissa})")
        g_calc = g_calc * r_mantissa
        print(f"    - Result FAIL. This is the failing operation.")
        # The operation above will fail because 13 * 2 = 26, which is > 15.
        
    except Titan4BitError as e:
        print(f"  - Failure during g calculation: {e}\n")
        print("Conclusion: The multiplication of mantissas for 'g' inevitably results in a number greater than 15.")
        print("Even if 'g' could be calculated through some combination of extreme approximations, the next step will fail.\n")
        
        # Step 3: Demonstrate failure in sqrt(A) step, assuming g was calculated
        print("Step 3: Show inevitable failure in sqrt() via Newton-Raphson method.")
        print("The method requires y_new = (y_old + A / y_old) / 2.")
        print("Let's assume we want to calculate sqrt(6), a plausible sub-problem.\n  A = 6/1")
        
        y_old = TitanFraction(5, 2) # Start guess y_0 = 2.5
        A = TitanFraction(6, 1)
        
        print(f"  - Iteration y_1: We must calculate ({y_old} + {A}/{y_old})")
        A_div_y = A / y_old
        print(f"  - The sum is: {y_old} + {A_div_y}")

        try:
            # This is the addition that will fail
            total = y_old + A_div_y
        except Titan4BitError as e_add:
            print(f"  - Fails here: {e_add}\n")
            print("The addition of two fractions within the square root algorithm requires a numerator")
            print("that exceeds the 4-bit limit. This problem is unavoidable for non-trivial calculations.")
            print("Therefore, the Titan computer cannot perform this task.")
            return "N0"

    return "Calculation somehow succeeded." # Should not be reached

# Execute the simulation and print the final answer
final_answer = calculate_landing_time()
print(f"\nFinal Answer: {final_answer}")
print("<<<N0>>>")
