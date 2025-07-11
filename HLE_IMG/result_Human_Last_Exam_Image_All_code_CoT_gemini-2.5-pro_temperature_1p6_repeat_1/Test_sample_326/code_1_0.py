import numpy as np

def solve_soliton_problem():
    """
    Solves the 2D Surface Flat-Top Soliton problem.
    """
    # Step 1: Define Material Parameters
    c_units = 1e11  # N/m^2
    d_units = 1e-11  # C/N
    eps0 = 9e-12    # F/m

    # KNbO3 (Orthorhombic, mm2)
    c_K = {'13': 0, '23': 0, '33': 3.0}
    d_K = {'31': -2.0, '32': 2.0, '33': 3.0}
    k_K = {'33': 300.0}
    rho_K = 4500.0

    # PMN-30%PT (Tetragonal, 4mm)
    c_P = {'13': 1.0, '33': 10.0/9.0}
    d_P = {'31': -100.0, '33': 200.0}
    k_P = {'33': 1000.0}
    rho_P = 8000.0

    # LiNbO3 (Trigonal, 3m)
    c_L = {'13': 3.0/4.0, '33': 5.0/2.0}
    d_L = {'31': -1.0/10.0, '33': 3.0/2.0}
    k_L = {'33': 30.0}
    rho_L = 4500.0

    # Step 2: Calculate longitudinal sound speed along the z-axis for each material.
    # This direction is assumed to correspond to the maximum velocity.

    # KNbO3: e_33 = d_31*c_13 + d_32*c_23 + d_33*c_33
    e33_K = (d_K['31']*c_K['13'] + d_K['32']*c_K['23'] + d_K['33']*c_K['33']) * c_units * d_units
    eps33_K = k_K['33'] * eps0
    c_prime_33_K = c_K['33']*c_units + e33_K**2 / eps33_K
    v_K = np.sqrt(c_prime_33_K / rho_K)

    # PMN-PT: Using tetragonal symmetry (d31=d32, c13=c23), e_33 = 2*d_31*c_13 + d_33*c_33
    e33_P = (2*d_P['31']*c_P['13'] + d_P['33']*c_P['33']) * c_units * d_units
    eps33_P = k_P['33'] * eps0
    c_prime_33_P = c_P['33']*c_units + e33_P**2 / eps33_P
    v_P = np.sqrt(c_prime_33_P / rho_P)

    # LiNbO3: Using trigonal symmetry (d31=d32, c13=c23), e_33 = 2*d_31*c_13 + d_33*c_33
    e33_L = (2*d_L['31']*c_L['13'] + d_L['33']*c_L['33']) * c_units * d_units
    eps33_L = k_L['33'] * eps0
    c_prime_33_L = c_L['33']*c_units + e33_L**2 / eps33_L
    v_L = np.sqrt(c_prime_33_L / rho_L)
    
    # Step 3: Determine v1 and v2.
    # The max longitudinal speed is max(v_K, v_P, v_L).
    # Since longitudinal modes are the fastest, v1 = v2 = max speed.
    v_max = max(v_K, v_P, v_L)
    v1 = v_max
    v2 = v_max

    # Step 4: Simplify NLSE and find M = max|Φ|
    # With v1=v2, the NLSE simplifies so that max|Φ|⁸ = ω.
    # The given logarithmic derivative implies ω = 17/324.
    omega_num = 17
    omega_den = 324
    omega = omega_num / omega_den

    max_phi = omega**(1.0/8.0)

    # Step 5: Calculate the final result 1 - max|Φ|.
    result = 1.0 - max_phi

    # Final output, showing the numbers in the equation for the answer.
    val1 = 1
    val2 = omega_num
    val3 = omega_den
    val4 = 8
    print(f"The calculation for the final answer is:")
    print(f"{val1} - ({val2}/{val3})**(1/{val4}) = {result:.12f}")

solve_soliton_problem()