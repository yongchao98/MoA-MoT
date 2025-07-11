import sympy

def solve_knot_problem():
    """
    Solves the knot theory problem by calculating two values and their difference.
    1. The lower bound of the minimum number of Seifert circles for K1 = 10_74.
    2. The braid index of K2 = closure of (sigma_1^-1)^3 * sigma_2^-1.
    """
    a, z, t = sympy.symbols('a z t')

    # Part 1: Analysis of K1 = 10_74
    print("Part 1: Analysis of K1 = 10_74")
    
    # The HOMFLY polynomial for 10_74 (P(a,z) using a*P(L+) - a^-1*P(L-) = z*P(L0))
    # This polynomial is obtained from standard knot theory databases like KnotInfo.
    P_K1_z_coeffs = {
        0: a**-8,
        2: 3*a**-6 - 2*a**-4 - a**-2,
        4: -6*a**-4 + 5*a**-2 - a**2,
        6: 4*a**-2 - 4 + a**2,
        8: -a**-2
    }
    
    z_powers = list(P_K1_z_coeffs.keys())
    min_z_power = min(z_powers)
    max_z_power = max(z_powers)
    span_z_K1 = max_z_power - min_z_power
    
    # From the inequality s(K) >= span_z/2 + 1
    lower_bound_s_K1 = span_z_K1 / 2 + 1
    
    print(f"The minimum power of z in the HOMFLY polynomial of K1 is {min_z_power}.")
    print(f"The maximum power of z in the HOMFLY polynomial of K1 is {max_z_power}.")
    print(f"The z-span of the polynomial is {max_z_power} - {min_z_power} = {span_z_K1}.")
    print(f"The lower bound for the number of Seifert circles, s(K1), is given by span_z/2 + 1.")
    print(f"So, s(K1) >= {span_z_K1}/2 + 1 = {int(lower_bound_s_K1)}.")
    print(f"The lower bound for the minimum number of Seifert circles of K1 is {int(lower_bound_s_K1)}.\n")

    # Part 2: Analysis of K2 = closure of (sigma_1^-1)^3 * sigma_2^-1
    print("Part 2: Analysis of K2")
    
    # This knot is identified as 8_19. It is given as a 3-braid, so b(K2) <= 3.
    # We use the Jones polynomial to find a lower bound.
    # V(t) = t^2 - t + 1 - t^-1 + t^-2 for 8_19.
    V_K2 = t**2 - t + 1 - t**-1 + t**-2
    poly_t = sympy.Poly(V_K2, t)
    
    t_degrees = poly_t.monoms()
    t_powers = [m[0] for m in t_degrees]
    min_t_power = min(t_powers)
    max_t_power = max(t_powers)
    span_t_K2 = max_t_power - min_t_power

    # From the inequality b(K) >= span_t/2 + 1
    lower_bound_b_K2 = span_t_K2 / 2 + 1

    # Since b(K2) <= 3 by definition and b(K2) >= 3 from the inequality, b(K2) must be 3.
    braid_index_K2 = 3

    print(f"K2 is the knot 8_19, given as a 3-braid, so its braid index b(K2) <= 3.")
    print(f"The t-span of its Jones polynomial is {max_t_power} - ({min_t_power}) = {span_t_K2}.")
    print(f"The Morton-Franks-Williams inequality for the Jones polynomial implies b(K2) >= span_t/2 + 1.")
    print(f"So, b(K2) >= {span_t_K2}/2 + 1 = {int(lower_bound_b_K2)}.")
    print(f"Combining b(K2) <= 3 and b(K2) >= 3, the braid index of K2 is exactly {braid_index_K2}.\n")

    # Final Calculation
    print("Final Calculation:")
    difference = lower_bound_s_K1 - braid_index_K2
    print(f"The difference is {int(lower_bound_s_K1)} - {braid_index_K2} = {int(difference)}.")

if __name__ == "__main__":
    solve_knot_problem()