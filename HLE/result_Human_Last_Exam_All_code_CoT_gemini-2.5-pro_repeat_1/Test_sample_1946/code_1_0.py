import numpy as np
from scipy.integrate import quad

def solve_toilet_paper_problem():
    """
    Calculates the time it takes for a falling toilet paper roll to unravel completely
    by applying the principle of conservation of energy and using numerical integration.
    """
    # --- Step 1: Define physical parameters and convert to SI units ---
    # Parameters from the problem description
    d_cyl_cm = 4.0
    t_paper_mm = 0.5
    N_wraps = 100
    m_paper_g = 200.0
    m_cyl_g = 20.0
    g = 9.8 # m/s^2

    # Convert to SI units (meters, kilograms)
    R_CYL = (d_cyl_cm / 100) / 2 # m
    T_PAPER = t_paper_mm / 1000 # m
    M_PAPER_TOTAL = m_paper_g / 1000 # kg
    M_CYL = m_cyl_g / 1000 # kg

    # --- Step 2: Calculate derived properties of the roll ---
    # The initial outer radius of the full roll
    R_OUTER_INITIAL = R_CYL + N_wraps * T_PAPER
    # The total length of the paper, calculated from the cross-sectional area
    # L_total = Area / thickness
    L_TOTAL = np.pi * (R_OUTER_INITIAL**2 - R_CYL**2) / T_PAPER
    # The linear mass density of the paper (mass per unit length)
    LAMBDA_PAPER = M_PAPER_TOTAL / L_TOTAL

    # --- Helper Function to calculate state-dependent properties as a function of y ---
    def get_state_properties(y):
        # y is the length of paper that has unrolled
        if y >= L_TOTAL: y = L_TOTAL
        
        # Remaining length of paper on the roll
        l_rem = L_TOTAL - y
        
        # Current outer radius of the roll
        r_outer_sq = R_CYL**2 + T_PAPER * l_rem / np.pi
        r_outer = np.sqrt(max(0, r_outer_sq))
        
        # Mass of the remaining paper on the roll
        m_paper_rem = LAMBDA_PAPER * l_rem
        
        # Total mass of the falling object (roll + remaining paper)
        M_total = M_CYL + m_paper_rem
        
        # Moment of inertia of the cardboard cylinder (as a thin shell)
        I_cyl = M_CYL * R_CYL**2
        # Moment of inertia of the remaining paper (as a hollow cylinder)
        I_paper_rem = 0.5 * m_paper_rem * (r_outer_sq + R_CYL**2)
        # Total moment of inertia
        I_total = I_cyl + I_paper_rem
        
        return r_outer, M_total, I_total

    # --- Step 3 & 4: Formulate the energy expressions and solve for velocity ---
    def velocity_sq(y):
        # v^2(y) is derived from KE_total(y) = -PE_total(y)
        if y <= 1e-12: return 0.0

        r, M, I = get_state_properties(y)
        if r <= 0: return 0.0

        # Numerator from potential energy: -PE = g*y*(M(y) + 0.5*lambda*y)
        numerator_v_sq = g * y * (M + 0.5 * LAMBDA_PAPER * y)
        
        # Denominator from kinetic energy factor: KE = v^2 * [ ... ]
        # where KE_factor = 0.5*(M(y) + I(y)/r^2) + (1/6)*lambda*y
        denominator_v_sq = 0.5 * (M + I / r**2) + (1/6.) * LAMBDA_PAPER * y
        
        if denominator_v_sq <= 0: return 0.0
        return numerator_v_sq / denominator_v_sq

    # --- Step 5: Define the integrand for calculating total time ---
    # Time t = integral(dt) = integral(dy/v)
    def integrand(y):
        v_sq = velocity_sq(y)
        if v_sq <= 0: return np.inf
        return 1.0 / np.sqrt(v_sq)

    # --- Step 6: Perform numerical integration ---
    # We integrate 1/v(y) from y=0 to y=L_total to find the time.
    time, _ = quad(integrand, 0, L_TOTAL, limit=200)

    # --- Print the explanation and result ---
    print("Step-by-step solution for the falling toilet paper problem:")
    print("\n1. Physical Parameters (SI Units):")
    print(f"   Gravity (g): {g:.1f} m/s^2")
    print(f"   Cylinder radius (r_c): {R_CYL:.3f} m")
    print(f"   Paper thickness (t): {T_PAPER:.4f} m")
    print(f"   Total paper mass (m_p): {M_PAPER_TOTAL:.3f} kg")
    print(f"   Cylinder mass (m_c): {M_CYL:.3f} kg")
    print(f"   Number of wraps (N): {N_wraps}")

    print("\n2. Derived Properties:")
    print(f"   Initial outer radius (R_0): {R_OUTER_INITIAL:.3f} m")
    print(f"   Total paper length (L): {L_TOTAL:.2f} m")
    print(f"   Paper linear density (lambda): {LAMBDA_PAPER:.4f} kg/m")

    print("\n3. Final Equation from Conservation of Energy:")
    print("   The square of the velocity (v) after a length 'y' has unrolled is given by:")
    print("   v(y)^2 = [g*y*(M(y) + 0.5*lambda*y)] / [0.5*(M(y) + I(y)/r(y)^2) + (1/6)*lambda*y]")
    print("   where M(y), I(y), and r(y) are the mass, moment of inertia, and radius of the roll.")
    
    print("\n4. Calculation:")
    print("   The total time is found by numerically integrating dt = dy/v(y) from y=0 to L.")

    print(f"\nFinal Result:")
    print(f"The time it takes for the toilet paper to reach the end of its roll is {time:.2f} seconds.")
    
    # Final answer in the required format
    print(f"\n<<<{time:.2f}>>>")

solve_toilet_paper_problem()