import math

def get_titan_representation(value_str):
    """A helper to format the Titan number representation for printing."""
    return value_str.replace("e", " * 10^")

def main():
    """
    Calculates the required rocket force using Titan 5-bit architecture rules.
    """
    print("[Titan Computer] Calculating required landing force F_rocket = F_g + F_d")
    print("-" * 60)

    # --- Step 1: Calculate Gravitational Force (F_g) ---
    # F_g = G * M * m / r^2
    # We use the most accurate 5-bit fractional approximations for constants.
    # G = 6.674e-11 => 20/3 * 10^-11
    # M_pandora is dominated by the shell. M ~ (4/3)*pi*rho_s*a_s^2*c_s
    # M ~ (4/3)*(22/7)*(3*10^2)*(2*10^6)^2*(2*10^6) ~ 1.0e22 kg => 10/1 * 10^21
    # m_probe = 50kg => 5/1 * 10^1
    # r = radius + altitude. Altitude is negligible. r ~ 2e6m => 2/1 * 10^6 => r^2 = 4/1 * 10^12
    
    print("[Analysis] Calculating gravitational force F_g = (G * M * m) / r^2")
    
    Fg_numerator_val = "(20/3 * 10^-11) * (10/1 * 10^21) * (5/1 * 10^1)"
    Fg_numerator_calc = "(20 * 10 * 5) / 3 * 10^(-11 + 21 + 1) = 1000/3 * 10^11"
    print(f"  Numerator (G*M*m) = {Fg_numerator_val}")
    print(f"  ... simplifies to = {get_titan_representation(Fg_numerator_calc)}")
    
    Fg_denominator_val = "(4/1 * 10^12)"
    print(f"  Denominator (r^2) = {Fg_denominator_val}")
    
    # F_g = (1000/3 * 10^11) / (4 * 10^12) = 1000/12 * 10^-1 = 250/3 * 10^-1 = 8.333... N
    # To represent 8.333..., the best Titan fraction is 25/3.
    Fg_final_num = 25
    Fg_final_den = 3
    Fg_final_exp = 0
    Fg_val = (Fg_final_num / Fg_final_den) * (10**Fg_final_exp)

    print(f"  Resulting F_g = (1000/12) * 10^-1 = 83.333 * 10^-1 = {Fg_val:.3f} N")
    print(f"  Final Titan representation for F_g = ({Fg_final_num} / {Fg_final_den}) * 10^{Fg_final_exp}")
    print("-" * 60)

    # --- Step 2: Calculate Deceleration Force (F_d) ---
    # F_d = (m * v^2) / (2*d)
    # m = 50 => 5/1 * 10^1
    # v = 300 => 3/1 * 10^2 => v^2 = 9/1 * 10^4
    # d = 5000 => 5/1 * 10^3 => 2d = 10/1 * 10^3
    print("[Analysis] Calculating deceleration force F_d = (m * v^2) / (2*d)")

    Fd_numerator_val = "(5/1 * 10^1) * (9/1 * 10^4)"
    Fd_numerator_calc = "45 * 10^5"
    print(f"  Numerator (m*v^2) = {Fd_numerator_val} = {get_titan_representation(Fd_numerator_calc)}")

    Fd_denominator_val = "(10/1 * 10^3)"
    print(f"  Denominator (2*d) = {Fd_denominator_val}")

    # F_d = (45 * 10^5) / (10 * 10^3) = 4.5 * 10^2 = 450 N
    # To represent 450, we use (450 / 10^3) * 10^3 = 0.45 * 10^3.
    # The fraction for 0.45 is 9/20.
    Fd_final_num = 9
    Fd_final_den = 20
    Fd_final_exp = 3
    Fd_val = (Fd_final_num / Fd_final_den) * (10**Fd_final_exp)

    print(f"  Resulting F_d = 4.5 * 10^2 = {Fd_val} N")
    print(f"  Final Titan representation for F_d = ({Fd_final_num} / {Fd_final_den}) * 10^{Fd_final_exp}")
    print("-" * 60)

    # --- Step 3: Final Summation ---
    # F_rocket = F_g + F_d = 8.333... + 450 = 458.333... N
    # We need to find the best Titan representation for 458.333...
    # Try with exp=3: 458.333... / 10^3 = 0.458333...
    # Best fraction for 0.458333... is 11/24.
    F_rocket_num = 11
    F_rocket_den = 24
    F_rocket_exp = 3
    F_rocket_val = (F_rocket_num / F_rocket_den) * (10**F_rocket_exp)

    print("[Analysis] Final summation F_rocket = F_g + F_d")
    print(f"  F_rocket = {Fg_val:.3f} N + {Fd_val} N = {Fg_val + Fd_val:.3f} N")
    print("Final Titan representation for F_rocket:")
    print(f"  ({Fg_final_num}/{Fg_final_den}) * 10^{Fg_final_exp} + ({Fd_final_num}/{Fd_final_den}) * 10^{Fd_final_exp} = ({F_rocket_num}/{F_rocket_den}) * 10^{F_rocket_exp}")
    print(f"  Final calculated force: {F_rocket_val:.3f} N")
    print("-" * 60)

    # --- Step 4: Error Calculation ---
    # High-precision calculation
    G_real = 6.674e-11
    rho_core = 1200
    r_core = 1e5
    M_core = (4/3) * math.pi * r_core**3 * rho_core
    rho_shell = 300
    a_shell = 2e6
    c_shell = 1.985e6
    M_shell = (4/3) * math.pi * (a_shell**2 * c_shell - r_core**3) * rho_shell
    M_pandora_real = M_core + M_shell
    m_probe_real = 50
    r_real = a_shell + 5000
    Fg_real = (G_real * M_pandora_real * m_probe_real) / (r_real**2)
    Fd_real = (50 * 300**2) / (2 * 5000)
    F_rocket_real = Fg_real + Fd_real
    
    error = abs(F_rocket_val - F_rocket_real)
    
    print(f"[Verification] Comparing against high-precision value")
    print(f"  Titan calculated value: {F_rocket_val:.3f} N")
    print(f"  High-precision value: {F_rocket_real:.3f} N")
    print(f"  Absolute Error = |{F_rocket_val:.3f} - {F_rocket_real:.3f}| = {error:.3f} N")
    
if __name__ == "__main__":
    main()
