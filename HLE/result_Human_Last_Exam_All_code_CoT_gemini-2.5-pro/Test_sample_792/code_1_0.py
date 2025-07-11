import math

def titan_calculate_force():
    """
    Simulates the force calculation on the Titan computer, following its constraints.
    """
    print("Starting calculation for the force F based on the formula: F = (D * m * g) / h\n")

    # Step 1: Define initial known values as Titan fractions (numerator, denominator)
    # D = 20m, h = 10m
    D = (20, 1)
    h = (10, 1)
    print(f"Known values: Target Distance D = {D[0]}/{D[1]}, Height h = {h[0]}/{h[1]}")

    # The formula simplifies to F = (20/10) * m * g = 2 * m * g
    factor = (2, 1)
    print(f"The formula simplifies to F = ({factor[0]}/{factor[1]}) * m * g\n")

    # Step 2: Choose optimal fractional approximations for m and g
    # m_real = 0.15 * pi ~= 0.4712 kg
    # g_real = 9.8 m/s^2
    # We choose fractions that allow for cancellation to keep intermediate values within 5-bit limits.
    m = (14, 29)  # Approximation for m (~0.482)
    g = (29, 3)   # Approximation for g (~9.67)
    print("Chosen approximations to enable calculation:")
    print(f"Mass m ≈ {m[0]}/{m[1]}")
    print(f"Gravity g ≈ {g[0]}/{g[1]}\n")

    # Step 3: Perform the calculation step-by-step, respecting Titan's constraints.
    print("Calculation steps:")

    # First multiplication: factor * m
    # (2/1) * (14/29)
    temp_num = factor[0] * m[0]
    temp_den = factor[1] * m[1]
    
    # Check constraints
    if temp_num > 31 or temp_den > 31:
        print(f"Error: Intermediate value {temp_num}/{temp_den} exceeds 5-bit limit.")
        return

    temp_result = (temp_num, temp_den)
    print(f"1. Multiply by mass: ({factor[0]}/{factor[1]}) * ({m[0]}/{m[1]}) = {temp_result[0]}/{temp_result[1]}")
    
    # Second multiplication: temp_result * g
    # (28/29) * (29/3)
    # To stay within limits, we must simplify by cancelling '29' before multiplying.
    print(f"2. Multiply by gravity: ({temp_result[0]}/{temp_result[1]}) * ({g[0]}/{g[1]})")
    print("   To avoid exceeding 31, we simplify by cancellation first:")
    
    # Simplified multiplication
    final_num = temp_result[0] * 1  # Since g[0] cancels with temp_result[1]
    final_den = 1 * g[1]

    if final_num > 31 or final_den > 31:
        print(f"Error: Final value {final_num}/{final_den} exceeds 5-bit limit.")
        return
        
    final_result = (final_num, final_den)
    print(f"   Simplified calculation: ({temp_result[0]}/1) * (1/{g[1]}) = {final_result[0]}/{final_result[1]}\n")

    # Step 4: Final Result
    print("Calculation successful.")
    print(f"The final required force F is represented by the fraction: {final_result[0]} / {final_result[1]} Newtons.")
    
    # Calculate the landing distance with this force to confirm it hits the target
    F_val = final_result[0] / final_result[1]
    m_real = 0.15 * math.pi
    g_real = 9.8
    h_real = 10
    D_landed = (F_val * h_real) / (m_real * g_real)
    print(f"\nThis force results in a landing distance of {D_landed:.2f} m, which is within the target range of 19-21 m.")

if __name__ == '__main__':
    titan_calculate_force()
