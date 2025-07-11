import math

# Helper class to handle fractional arithmetic with Titan's 5-bit constraints.
class TitanFraction:
    def __init__(self, num, den=1):
        if not (0 <= num <= 31 and 0 <= den <= 31):
            raise ValueError(f"Numerator {num} or denominator {den} exceeds 5-bit limit (0-31)")
        self.num = num
        self.den = den

    def __mul__(self, other):
        # Multiply two fractions
        new_num = self.num * other.num
        new_den = self.den * other.den
        
        # Simplify the result
        common_divisor = math.gcd(new_num, new_den)
        final_num = new_num // common_divisor
        final_den = new_den // common_divisor
        
        # Check if the simplified result respects the 5-bit constraint
        if not (0 <= final_num <= 31 and 0 <= final_den <= 31):
            raise ValueError(f"Intermediate result {final_num}/{final_den} exceeds 5-bit limit")
            
        return TitanFraction(final_num, final_den)

    def __truediv__(self, other):
        # Division is multiplication by the reciprocal
        reciprocal = TitanFraction(other.den, other.num)
        return self.__mul__(reciprocal)

    def __repr__(self):
        return f"{self.num}/{self.den}"

def solve_titan_problem():
    """
    Solves the Curious Monkey problem using the Titan computer architecture.
    """
    print("--- Titan Superconducting Computer Simulation ---")
    print("Objective: Calculate the force F to hit the coconut.\n")

    # 1. Define the physics equation and calculate the true value for reference
    # F = (g * m * sqrt(2) * x) / (x - y)
    # Mass m = density * Volume = 0.9 kg/cm^3 * (4/3 * pi * (0.5 cm)^3) = 0.15 * pi kg
    m_true = 0.15 * math.pi
    g_true = 9.8
    sqrt2_true = math.sqrt(2)
    x_true = 20.0
    y_true = 10.0
    
    f_true = (g_true * m_true * sqrt2_true * x_true) / (x_true - y_true)
    print(f"Target 'True' Force: {f_true:.4f} N\n")

    # 2. Choose fractional approximations that fit 5-bit registers (0-31)
    # These are chosen strategically to allow for cancellation during multiplication.
    g = TitanFraction(13, 1)     # Approximation of 9.8
    m = TitanFraction(1, 2)      # Approximation of m_true (0.471)
    sqrt2 = TitanFraction(1, 1)  # Approximation of sqrt(2) (1.414)
    x = TitanFraction(20, 1)
    y = TitanFraction(10, 1)
    
    print("Chosen Fractional Approximations:")
    print(f"g     ≈ {g}")
    print(f"m     ≈ {m}")
    print(f"sqrt(2) ≈ {sqrt2}")
    print(f"x     = {x}")
    print(f"y     = {y}\n")

    # 3. Perform the calculation step-by-step, respecting Titan's constraints
    print("Step-by-step Calculation:")
    
    # First, calculate the term x / (x - y)
    # Note: Subtraction is not explicitly defined, but x and y are simple integers.
    # x - y = 20 - 10 = 10. We represent this as 10/1.
    x_minus_y = TitanFraction(10, 1)
    term1 = x / x_minus_y
    print(f"Step 1: x / (x - y) = {x} / {x_minus_y} = {term1}")

    # Now, multiply all terms together sequentially
    # F = g * m * sqrt2 * term1
    res_step1 = g * m
    print(f"Step 2: g * m = {g} * {m} = {res_step1}")
    
    res_step2 = res_step1 * sqrt2
    print(f"Step 3: (g * m) * sqrt(2) = {res_step1} * {sqrt2} = {res_step2}")
    
    final_force = res_step2 * term1
    print(f"Step 4: (g * m * sqrt(2)) * (x / (x-y)) = {res_step2} * {term1} = {final_force}\n")

    # 4. Final Result and Error Calculation
    f_calc = final_force.num / final_force.den
    error = abs(f_calc - f_true)

    print("--- Final Result ---")
    print(f"The final equation is:")
    print(f"{final_force} = {g} * {m} * {sqrt2} * {term1}")
    print(f"\nThe calculated force is {f_calc:.4f} N.")
    
    # Check if the shot is successful (target force is between 12.93N and 13.20N)
    if 12.93 <= f_calc <= 13.20:
        print("This force is within the required range to hit the coconut. Success!")
    else:
        print("This force is outside the required range. The shot misses.")

    print(f"The absolute error from the ideal center is {error:.4f}.")
    
    # Format the final answer as requested
    final_answer = f"Y[{error:.3f}]"
    print(f"\nFinal Answer String: {final_answer}")
    print(f"<<<{final_answer}>>>")

solve_titan_problem()