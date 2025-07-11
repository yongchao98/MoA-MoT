import math

def g_hard_rods(x, eta):
    """
    Calculates the exact radial distribution function g(x) for 1D hard rods.
    x is the distance in units of sigma (r/sigma).
    eta is the packing fraction.
    """
    if x < 1:
        return 0.0
    
    lambda_val = eta / (1.0 - eta)
    n_max = math.floor(x)
    
    total_sum = 0.0
    for n in range(1, n_max + 1):
        try:
            term = (math.pow(lambda_val, n) / math.factorial(n - 1)) * \
                   math.pow(x - n, n - 1) * \
                   math.exp(-lambda_val * (x - n))
            total_sum += term
        except ValueError: # handles (x-n)^(n-1) where base is 0 and exp is < 0
            if x == n and n-1 == 0:
                 term = (math.pow(lambda_val, n) / math.factorial(n - 1)) * \
                   1 * \
                   math.exp(-lambda_val * (x - n))
                 total_sum += term
            else:
                 continue


    return total_sum / eta

def solve():
    # Step 1-4: Identification of plots based on physical characteristics
    # g(r) for {SS, SR, R, HS, TW}
    g_ss = 1  # Repulsive shoulder from r=1 to 1.5 creates a dip.
    g_sr = 9  # Sticky potential creates a delta-function-like peak at contact (r=1).
    g_r = 5   # Repulsive ramp potential causes g(r) to ramp up from r=1 to 1.5.
    g_hs = 3  # Hard spheres have features only at integer multiples of sigma. By elimination.
    g_tw = 7  # Attractive triangle well creates a large peak at contact (r=1).
    
    g_r_indices = [g_ss, g_sr, g_r, g_hs, g_tw]
    
    # S(k) for {SS, SR, R, HS, TW}
    # S(0) is related to compressibility. Attraction increases S(0), repulsion decreases it.
    # Order of S(0): SR > TW > 1 > HS > R > SS
    # From plots: S(0)_2=4 > S(0)_6~2.5 > 1 > S(0)_8~0.2 > S(0)_4~0
    # SR: S(0) = 1/(1-alpha)^2 = 1/(1-1.5)^2 = 4. This matches plot 2.
    s_sr = 2
    # TW: Attractive, so S(0) > 1. This matches plot 6.
    s_tw = 6
    # SS: Most repulsive potential (on average), so lowest S(0). Matches plot 4.
    s_ss = 4
    # R: Repulsive, S(0) < S(0)_HS. Matches plot 8.
    s_r = 8
    # HS: S(0)=(1-eta)^2 = (2/3)^2 ~ 0.44. No plot matches. Thus, HS is the unique system.
    s_hs = 0
    
    s_k_indices = [s_ss, s_sr, s_r, s_hs, s_tw]

    # Step 5: Calculate R_max for the unique system (HS)
    eta = 1.0/3.0
    r_values = [i / 2.0 for i in range(3, 22, 2)] # r/sigma = {1.5, 2.5, ..., 10.5}
    
    ratios = []
    print("Calculating R_max for the unique system (Hard Spheres, HS):")
    print("R_g(r/sigma) = g(r/sigma + 1) / g(r/sigma)")
    print("-" * 50)
    
    for r_over_sigma in r_values:
        g_r = g_hard_rods(r_over_sigma, eta)
        g_r_plus_1 = g_hard_rods(r_over_sigma + 1, eta)
        ratio = g_r_plus_1 / g_r
        ratios.append(ratio)
        print(f"r/sigma = {r_over_sigma:4.1f}: R_g({r_over_sigma:4.1f}) = g({r_over_sigma+1:4.1f})/g({r_over_sigma:4.1f}) = {g_r_plus_1:.6f} / {g_r:.6f} = {ratio:.6f}")

    r_max = max(ratios)
    
    # Step 6: Final Answer
    final_sequence = g_r_indices + s_k_indices + [r_max]
    
    print("-" * 50)
    print("The final sequence {g(SS), g(SR), g(R), g(HS), g(TW), S(SS), S(SR), S(R), S(HS), S(TW), R_max} is:")
    # Create the string representation for the final list of numbers.
    # We use a special format for R_max to show more precision.
    result_str = "{"
    for i in range(10):
        result_str += str(final_sequence[i]) + ", "
    result_str += f"{final_sequence[10]:.17f}" + "}"
    print(result_str)

solve()
<<< {1, 9, 5, 3, 7, 4, 2, 8, 0, 6, 1.00235339620894560}>>>