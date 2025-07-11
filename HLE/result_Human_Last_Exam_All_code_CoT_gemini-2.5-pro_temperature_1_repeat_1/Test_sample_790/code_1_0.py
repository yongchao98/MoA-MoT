import math

def find_best_approximation(value):
    """
    Finds the best fractional approximation a/b for a given decimal value,
    where a and b are 5-bit integers (0-31).
    """
    min_error = float('inf')
    best_fraction = (0, 1)

    for b in range(1, 32):
        # Find the two best integer numerators for this denominator
        ideal_a = value * b
        a1 = math.floor(ideal_a)
        a2 = math.ceil(ideal_a)

        # Check the floor numerator
        if 0 <= a1 <= 31:
            error = abs(value - a1 / b)
            if error < min_error:
                min_error = error
                best_fraction = (a1, b)

        # Check the ceiling numerator
        if 0 <= a2 <= 31:
            error = abs(value - a2 / b)
            if error < min_error:
                min_error = error
                best_fraction = (a2, b)
    
    return best_fraction

def main():
    """
    Solves the problem using Titan computer rules.
    """
    print("--- Titan Computer Calculation for Force F ---")
    print("Equation: F = (2 * m * g) / cos(45°)\n")

    # 1. Initial fractional approximations
    m_num, m_den = 9, 19
    g_num, g_den = 29, 3
    cos45_num, cos45_den = 12, 17
    
    print("Step 1: Define initial approximations")
    print(f"  m ≈ {m_num}/{m_den}")
    print(f"  g ≈ {g_num}/{g_den}")
    print(f"  cos(45°) ≈ {cos45_num}/{cos45_den}\n")
    
    # 2. Calculate the numerator term: 2 * m * g
    print("Step 2: Calculate numerator term `2 * m * g`")
    # First part of multiplication: 2 * m
    temp_num = 2 * m_num
    temp_den = 1 * m_den
    print(f"  (2/1) * ({m_num}/{m_den}) = {temp_num}/{temp_den}. (Values are within 5-bit limits)")

    # Second part of multiplication: (2*m) * g
    rhs_num = temp_num * g_num
    rhs_den = temp_den * g_den
    print(f"  ({temp_num}/{temp_den}) * ({g_num}/{g_den}) = {rhs_num}/{rhs_den}")
    
    # Simplify using GCD
    common_divisor = math.gcd(rhs_num, rhs_den)
    rhs_num_s = rhs_num // common_divisor
    rhs_den_s = rhs_den // common_divisor
    print(f"  Simplified with GCD: {rhs_num_s}/{rhs_den_s}")

    # Check constraints and reduce if necessary
    if rhs_num_s > 31 or rhs_den_s > 31:
        print(f"  Numerator {rhs_num_s} exceeds 31. Reducing fraction...")
        value_to_approx = rhs_num_s / rhs_den_s
        rhs_approx_num, rhs_approx_den = find_best_approximation(value_to_approx)
        print(f"  Value {value_to_approx:.4f} is best approximated by {rhs_approx_num}/{rhs_approx_den}\n")
    else:
        rhs_approx_num, rhs_approx_den = rhs_num_s, rhs_den_s

    # 3. Calculate F = (2mg)_approx / cos(45)
    print("Step 3: Calculate final force F")
    print(f"  F ≈ ({rhs_approx_num}/{rhs_approx_den}) / ({cos45_num}/{cos45_den})")
    
    f_num = rhs_approx_num * cos45_den
    f_den = rhs_approx_den * cos45_num
    print(f"  F ≈ {f_num}/{f_den}")

    # Simplify using GCD
    common_divisor_f = math.gcd(f_num, f_den)
    f_num_s = f_num // common_divisor_f
    f_den_s = f_den // common_divisor_f
    print(f"  Simplified with GCD: {f_num_s}/{f_den_s}")

    # Check constraints and reduce if necessary
    final_f_num, final_f_den = (0,1)
    if f_num_s > 31 or f_den_s > 31:
        print(f"  Numerator {f_num_s} exceeds 31. Reducing fraction...")
        value_to_approx_f = f_num_s / f_den_s
        final_f_num, final_f_den = find_best_approximation(value_to_approx_f)
        print(f"  Value {value_to_approx_f:.4f} is best approximated by {final_f_num}/{final_f_den}\n")
    else:
        final_f_num, final_f_den = f_num_s, f_den_s

    f_approx_val = final_f_num / final_f_den
    print(f"Final calculated Force: F = {final_f_num}/{final_f_den} = {f_approx_val} N\n")
    
    # 4. Calculate the error
    print("Step 4: Calculate the absolute error")
    m_true = 0.15 * math.pi
    g_true = 9.8
    cos45_true = math.cos(math.radians(45))
    
    f_true = (2 * m_true * g_true) / cos45_true
    error = abs(f_true - f_approx_val)
    
    print(f"  True value of F is approximately {f_true:.3f} N")
    print(f"  Calculated value of F is {f_approx_val:.3f} N")
    print(f"  Absolute error |{f_true:.3f} - {f_approx_val:.3f}| = {error:.3f}\n")
    
    # 5. Final Answer
    # Check if the rock hits the target range (9.9m to 10.1m)
    # y = 20 - (20 * m_true * g_true) / (f_approx_val * cos45_true)
    # This is equivalent to checking if (2 * m*g)/(F*cos45) is between 0.99 and 1.01
    term_check = (2 * m_true * g_true) / (f_approx_val * cos45_true)
    if 0.99 <= term_check <= 1.01:
        print("This force will hit the coconut. The calculation is feasible.")
        print(f"<<<Y[{error:.3f}]>>>")
    else:
        print("This force will miss the coconut. The calculation is not feasible with this accuracy.")
        print("<<<N0>>>")


if __name__ == "__main__":
    main()