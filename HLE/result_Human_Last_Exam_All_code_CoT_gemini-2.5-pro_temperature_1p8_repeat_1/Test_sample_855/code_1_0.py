import math

def find_best_fraction(value, max_val=31):
    """
    Finds the best fractional approximation n/d for a given value,
    where n and d are <= max_val.
    """
    best_n, best_d = 1, 1
    min_error = abs(value - 1.0)

    for d in range(1, max_val + 1):
        n = int(round(value * d))
        if n > max_val or n == 0:
            continue
        
        error = abs(value - n / d)
        if error < min_error:
            min_error = error
            best_n, best_d = n, d
            
    return best_n, best_d

def represent_as_titan(value):
    """
    Represents a floating-point value in the optimal Titan format (n/d) * 10^exp.
    """
    if value == 0:
        return 0, 1, 0
    
    # Choose an exponent so the mantissa is between 1 and 10
    exponent = int(math.floor(math.log10(abs(value))))
    mantissa = value / (10**exponent)
    
    # Sometimes, shifting the exponent can yield a better fraction.
    # We test the original exponent and one lower.
    mantissa1 = mantissa
    exp1 = exponent
    n1, d1 = find_best_fraction(mantissa1)
    error1 = abs(mantissa1 - n1/d1)

    mantissa2 = mantissa * 10
    exp2 = exponent - 1
    n2, d2 = find_best_fraction(mantissa2)
    error2 = abs(mantissa2 - n2/d2)

    # Return the representation with the smaller error
    if error1 < error2:
        return n1, d1, exp1
    else:
        return n2, d2, exp2

def main():
    """
    Performs the calculation for the Pioneer landing force on Pandora using
    the simulated Titan 5-bit architecture.
    """
    print("[Titan Calculation Steps]")
    
    # --- Part 1: Constants in Titan and Real Format ---
    # Titan uses fractional approximations.
    G_n, G_d, G_exp = 20, 3, -11  # Approx for 6.674e-11
    m_n, m_d, m_exp = 5, 1, 1      # For 50 kg
    v_n, v_d, v_exp = 3, 1, 2      # For 300 m/s
    d_n, d_d, d_exp = 5, 1, 3      # For 5000 m
    # For Pandora's Mass (M) and Radius (r), we pre-calculate and approximate
    # M = 9.6e21 kg -> 9.5 -> 19/2. So M is (19/2) * 10^21
    M_n, M_d, M_exp = 19, 2, 21
    # r = 2,005,000 m. Approx as 2e6. r^2 = 4e12.
    r2_n, r2_d, r2_exp = 4, 1, 12

    # --- Part 2: Calculate F_gravity = G * M * m / r^2 ---
    # Intermediate values are re-represented in Titan format at each step.
    # 1. G * M
    val_gm = (G_n/G_d * 10**G_exp) * (M_n/M_d * 10**M_exp)
    gm_n, gm_d, gm_exp = represent_as_titan(val_gm)

    # 2. (G * M) * m
    val_gmm = (gm_n/gm_d * 10**gm_exp) * (m_n/m_d * 10**m_exp)
    gmm_n, gmm_d, gmm_exp = represent_as_titan(val_gmm)

    # 3. F_gravity = (G * M * m) / r^2
    val_fg = (gmm_n/gmm_d * 10**gmm_exp) / (r2_n/r2_d * 10**r2_exp)
    fg_n, fg_d, fg_exp = represent_as_titan(val_fg)
    
    print("1. Gravitational Force (F_gravity):")
    print(f"   - Calculated as (G*M*m)/r^2")
    print(f"   - Titan Representation: ({fg_n} / {fg_d}) * 10^{fg_exp} N")
    
    # --- Part 3: Calculate F_decel = (m * v^2) / (2d) ---
    # 1. v^2
    val_v2 = ((v_n/v_d * 10**v_exp))**2
    v2_n, v2_d, v2_exp = represent_as_titan(val_v2)

    # 2. m * v^2
    val_mv2 = (m_n/m_d * 10**m_exp) * (v2_n/v2_d * 10**v2_exp)
    mv2_n, mv2_d, mv2_exp = represent_as_titan(val_mv2)

    # 3. 2 * d
    val_2d = 2 * (d_n/d_d * 10**d_exp)
    d2_n, d2_d, d2_exp = represent_as_titan(val_2d)

    # 4. F_decel = (m * v^2) / (2d)
    val_fd = (mv2_n/mv2_d * 10**mv2_exp) / (d2_n/d2_d * 10**d2_exp)
    fd_n, fd_d, fd_exp = represent_as_titan(val_fd)
    
    print("\n2. Deceleration Force (F_decel):")
    print(f"   - Calculated as (m*v^2)/(2d)")
    print(f"   - Titan Representation: ({fd_n} / {fd_d}) * 10^{fd_exp} N")
    
    # --- Part 4: Calculate Total Rocket Force ---
    # F_rocket = F_gravity + F_decel. Add values, then find best Titan representation.
    val_total = val_fg + val_fd
    total_n, total_d, total_exp = represent_as_titan(val_total)
    
    print("\n3. Total Rocket Force (F_rocket = F_gravity + F_decel):")
    print("   - F_gravity is much smaller than F_decel, so the sum is approximately F_decel.")
    print("   - A precise fractional sum is not possible within constraints.")
    print("   - We find the best Titan approximation for the combined numerical value.")
    print("\nFinal Equation:")
    print(f"F_rocket = ({fg_n}/{fg_d})*10^{fg_exp} + ({fd_n}/{fd_d})*10^{fd_exp} "
          f"≈ ({total_n} / {total_d}) * 10^{total_exp} N")
    titan_result = (total_n/total_d) * 10**total_exp
    print(f"           = {val_fg:.3f} N + {val_fd:.1f} N = {val_total:.3f} N")
    print(f"           Approximated as: {titan_result:.1f} N")

    # --- Part 5: Calculate Error ---
    G_real = 6.674e-11
    M_real = (1.2e3 * 4/3 * math.pi * (1e5)**3) + \
             (0.3e3 * (4/3 * math.pi * (2e6)**2 * 1.985e6 - 4/3 * math.pi * (1e5)**3))
    m_real = 50.0
    r_real = 2.005e6 # At 5000m altitude
    d_real = 5000.0
    v_real = 300.0

    fg_real = G_real * M_real * m_real / (r_real**2)
    fd_real = (m_real * v_real**2) / (2 * d_real)
    f_total_real = fg_real + fd_real
    
    absolute_error = abs(titan_result - f_total_real)
    
    print(f"\nReal value (using floats): {f_total_real:.3f} N")
    print(f"Absolute error: |{titan_result:.1f} - {f_total_real:.3f}| = {absolute_error:.3f} N")
    print(f"\nCan you use Titan to calculate this force? Yes.")
    print(f"The smallest absolute error we can produce is {absolute_error:.3f} N.")
    
    # Final answer format for the system
    final_answer = f"Y[{absolute_error:.3f}]"
    print("\n<<<" + final_answer + ">>>")

if __name__ == '__main__':
    main()
