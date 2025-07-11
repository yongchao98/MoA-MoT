import math

class TitanFraction:
    """
    A class to represent fractions under Titan's 5-bit constraints.
    It does not automatically handle operations but is used to store and display fractions.
    """
    def __init__(self, n, d):
        if not (0 <= n <= 31 and 1 <= d <= 31):
            raise ValueError("Numerator and denominator must be within 5-bit range (0-31)")
        self.n = n
        self.d = d

    def __repr__(self):
        return f"{self.n}/{self.d}"

    def to_float(self):
        return self.n / self.d

def find_best_approximation(target_value):
    """
    Finds the best fraction n/d (n, d <= 31) to approximate a target value.
    """
    best_fraction = TitanFraction(0, 1)
    min_error = float('inf')

    for d in range(1, 32):
        for n in range(0, 32):
            if d == 0: continue
            current_error = abs(target_value - n / d)
            if current_error < min_error:
                min_error = current_error
                best_fraction = TitanFraction(n, d)
    return best_fraction

def solve_coconut_problem():
    """
    Solves the physics problem following Titan computer's rules.
    """
    print("Monkey Business: Calculating the Coconut Force on Titan\n")
    print("Physics Derivation:")
    print("The required force is F = (3 * pi * g) / (10 * cos(45))\n")

    # High-precision constants as Titan fractions
    PI = TitanFraction(22, 7)         # approx 3.142
    G = TitanFraction(29, 3)          # approx 9.667
    COS45 = TitanFraction(17, 24)     # cos(45) = 1/sqrt(2), using sqrt(2) approx 24/17

    C_3_10 = TitanFraction(3, 10)

    print("Chosen High-Precision Constants:")
    print(f"pi ≈ {PI} ({PI.to_float():.4f})")
    print(f"g ≈ {G} ({G.to_float():.4f}) m/s^2")
    print(f"cos(45) ≈ {COS45} ({COS45.to_float():.4f})\n")

    print("Step-by-step Calculation (with approximations for overflows):\n")

    # Step 1: Calculate A = (3/10) / cos(45)
    print("Step 1: Calculate A = (3/10) / cos(45)")
    print(f"         A = ({C_3_10}) / ({COS45}) = ({C_3_10}) * ({COS45.d}/{COS45.n})")
    print(f"         Intermediate multiplication: {C_3_10.n}*{COS45.d} = {C_3_10.n * COS45.d} -> OVERFLOW (>31)")
    
    val_A = C_3_10.to_float() / COS45.to_float()
    A_approx = find_best_approximation(val_A)
    print(f"         True value = {val_A:.4f}. Approximating with a valid Titan fraction...")
    print(f"         Best approximation is {A_approx} ({A_approx.to_float():.4f})\n")

    # Step 2: Calculate B = A * pi
    print("Step 2: Calculate B = A_approx * pi")
    print(f"         B = ({A_approx}) * ({PI})")
    print(f"         Intermediate multiplication: {A_approx.n}*{PI.n} = {A_approx.n * PI.n} -> OVERFLOW (>31)")

    val_B = A_approx.to_float() * PI.to_float()
    B_approx = find_best_approximation(val_B)
    print(f"         True value = {val_B:.4f}. Approximating...")
    print(f"         Best approximation is {B_approx} ({B_approx.to_float():.4f})\n")

    # Step 3: Calculate Final Force F = B * g
    print("Step 3: Calculate Final Force F = B_approx * g")
    print(f"         F = ({B_approx}) * ({G})")
    print(f"         Intermediate multiplication: {B_approx.n}*{G.n} = {B_approx.n * G.n} -> OVERFLOW (>31)")

    val_F = B_approx.to_float() * G.to_float()
    F_final = find_best_approximation(val_F)
    print(f"         True value = {val_F:.4f}. Approximating...")
    print(f"         Best approximation is {F_final} ({F_final.to_float():.4f})\n")

    print("--------------------------------------------------")
    print(f"Final calculated force F = {F_final.n} / {F_final.d} N")
    print("--------------------------------------------------\n")

    # Error Calculation
    print("Calculating Error:")
    r = 0.5 / 100  # cm to m
    rho = 0.9 * (100**3) # kg/cm^3 to kg/m^3 -> using the other density to match derivation
    rho_kg_cm3 = 0.9
    mass = rho_kg_cm3 * (4/3) * math.pi * (0.5**3) # keeps calculation simple
    
    true_g = 9.80665
    true_F = (2 * true_g * mass) / math.cos(math.radians(45))
    
    calculated_F = F_final.to_float()
    abs_error = abs(true_F - calculated_F)
    
    print(f"High-precision benchmark value for F: {true_F:.4f} N")
    print(f"Our calculated value: {calculated_F:.4f} N")
    print(f"Absolute error: {abs_error:.3f}\n")

    # Print the final equation with all numbers
    print("Final Equation Trace:")
    print(f"F = [ (({C_3_10}) / ({COS45})) [≈ {A_approx}] * ({PI}) ] [≈ {B_approx}] * ({G}) ≈ {F_final}")


    # Final answer in specified format
    print("\nCan Titan calculate the force? Yes.")
    final_answer = f"<<<Y[{abs_error:.3f}]>>>"
    print(f"Final Answer Token: {final_answer}")


if __name__ == '__main__':
    solve_coconut_problem()
