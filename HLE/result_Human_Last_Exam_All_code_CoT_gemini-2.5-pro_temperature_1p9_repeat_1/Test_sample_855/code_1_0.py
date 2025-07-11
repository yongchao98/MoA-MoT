import math

class TitanNumber:
    """
    A helper class to represent numbers in Titan's format.
    It stores a number as a fraction (n/d) and a base-10 exponent.
    n and d are constrained to be 5-bit integers (0-31).
    """
    def __init__(self, n, d, exp=0):
        # In a full simulation, checks would prevent n, d from exceeding 31.
        self.n = n
        self.d = d
        self.exp = exp

    def value(self):
        """Calculates the decimal value of the Titan number."""
        if self.d == 0:
            return float('inf')
        return (self.n / self.d) * (10 ** self.exp)

    def __str__(self):
        """String representation of the Titan number."""
        return f"({self.n}/{self.d}) * 10^{self.exp}"

def find_best_titan_approximation(target_value):
    """
    Finds the best TitanNumber representation for a given decimal value.
    This simulates the processor finding the most accurate way to store a calculation result.
    """
    best_approx = None
    min_err = float('inf')

    # We check a reasonable range of exponents for a force in Newtons.
    for exp_candidate in range(-2, 4):
        mantissa = target_value / (10**exp_candidate)
        
        # For the resulting mantissa, find the best fraction n/d where n,d <= 31.
        for d_candidate in range(1, 32):
            n_candidate = round(mantissa * d_candidate)
            if 0 <= n_candidate <= 31:
                current_num = TitanNumber(n_candidate, d_candidate, exp_candidate)
                err = abs(current_num.value() - target_value)
                if err < min_err:
                    min_err = err
                    best_approx = current_num

    return best_approx

def solve_landing_problem():
    """
    Main function to solve the problem and output the reasoning and final answer.
    """
    # --- Step 1: Precise 'Real-World' Calculation (Ground Truth) ---
    G_real = 6.67430e-11
    m_probe = 50.0
    v_i = 300.0
    landing_dist = 5000.0
    r_core_m, rho_core_kg_m3 = 100e3, 1200.0
    r_eq_m, r_pol_m = 2000e3, 1985e3
    rho_shell_kg_m3 = 300.0

    V_core = (4/3) * math.pi * r_core_m**3
    M_core = V_core * rho_core_kg_m3
    V_ellipsoid = (4/3) * math.pi * r_eq_m**2 * r_pol_m
    V_shell = V_ellipsoid - V_core
    M_shell = V_shell * rho_shell_kg_m3
    M_pandora = M_core + M_shell
    
    r_probe_center = r_eq_m + landing_dist
    F_g_real = (G_real * M_pandora * m_probe) / (r_probe_center**2)
    acceleration_needed = v_i**2 / (2 * landing_dist)
    F_d_real = m_probe * acceleration_needed
    F_total_real = F_d_real + F_g_real

    # --- Step 2: Titan Architecture Calculation ---
    # We must calculate F_rocket = F_d + F_g under Titan's rules.
    
    # 2a: Deceleration Force (F_d)
    # F_d = m * a = 50 * (300^2 / (2*5000)) = 450 N.
    # We can represent 450 as 4.5 * 10^2. 4.5 can be the fraction 9/2.
    F_d_titan = TitanNumber(9, 2, 2)
    
    # 2b: Gravitational Force (F_g)
    # Direct calculation of F_g proves impossible as intermediate products (e.g., in calculating
    # mass) result in numerators or denominators larger than 31.
    # Therefore, we must use an approximation. The real gravitational acceleration g is ~0.1656 m/s^2.
    # The fraction 1/6 (~0.1667) is a very close and simple approximation.
    # So we approximate g_titan = 1/6.
    # F_g = m_probe * g = 50 * (1/6) = 25/3. This is a valid TitanNumber.
    F_g_titan = TitanNumber(25, 3, 0)
    
    # 2c: Total Force and Final Approximation
    # The sum is 450 + 25/3 = 1375/3 ≈ 458.333 N.
    # This result must be stored as a valid TitanNumber. We find the best approximation.
    titan_sum_value = F_d_titan.value() + F_g_titan.value()
    F_total_titan = find_best_titan_approximation(titan_sum_value)
    
    final_error = abs(F_total_titan.value() - F_total_real)
    
    # --- Step 3: Print Final Output ---
    print("The final equation is F_rocket = F_deceleration + F_gravity.")
    print("Here is the calculation using Titan's fractional arithmetic:\n")

    print(f"1. Deceleration force F_d = {F_d_titan.value()} N, represented as: {F_d_titan.n} / {F_d_titan.d} * 10^{F_d_titan.exp}")
    print(f"2. Gravitational force F_g ≈ {F_g_titan.value():.3f} N, represented as: {F_g_titan.n} / {F_g_titan.d} * 10^{F_g_titan.exp}")

    print("\n3. Adding these forces gives the total required rocket force:")
    # This is the key line demonstrating the final calculation as requested
    print(f"({F_d_titan.n} / {F_d_titan.d} * 10^{F_d_titan.exp}) + ({F_g_titan.n} / {F_g_titan.d}) = {F_total_titan.value()} N (approximated)")
    
    print(f"\nThe precise value is {F_total_real:.3f} N.")
    print(f"The Titan computer's best approximation is {F_total_titan.value()} N.")
    print(f"The smallest absolute error is |{F_total_titan.value()} - {F_total_real:.3f}| = {final_error:.3f} N.")

    answer = f"Y{final_error:.3f}"
    print(f"\n<<<Y{final_error:.3f}>>>")

solve_landing_problem()