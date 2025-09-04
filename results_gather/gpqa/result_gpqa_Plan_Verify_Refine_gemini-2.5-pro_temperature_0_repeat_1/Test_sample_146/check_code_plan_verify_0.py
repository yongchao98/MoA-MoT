import numpy as np

def check_annihilation_velocity():
    """
    This function checks the correctness of the provided answer for the particle annihilation problem.
    It recalculates the velocity of particle A from first principles and compares it to the given answer.
    """

    # --- Define Constants and Problem Parameters ---
    # Rest mass energy of a proton (and antiproton) in MeV.
    # Using a standard value from CODATA.
    m_p_c2 = 938.27208816  # MeV
    
    # Rest mass energy of particle A, as given in the question.
    m_A_c2 = 300.0  # MeV
    
    # The answer to check is Option B, which corresponds to v = 0.77c.
    provided_answer_v_over_c = 0.77

    # --- Step 1: Apply Conservation of Energy ---
    # The problem states the antiproton is "slowly moving", which implies we can
    # approximate the initial kinetic energy of the system as zero.
    # The initial energy (E_initial) is the sum of the rest mass energies of the proton and antiproton.
    E_initial = 2 * m_p_c2

    # The final state consists of four A particles (2 A+ and 2 A-).
    # To conserve momentum from an initial state of zero momentum, the final particles must
    # move in such a way that the total momentum is zero. By symmetry, all four particles
    # will have the same magnitude of velocity and thus the same energy.
    # The total energy of one A particle is E_A = gamma * m_A * c^2.
    # The final energy (E_final) is the sum of the energies of the four A particles.
    # E_final = 4 * gamma * m_A * c^2

    # According to the principle of conservation of energy, E_initial = E_final.
    # 2 * m_p_c2 = 4 * gamma * m_A_c2

    # --- Step 2: Solve for the Lorentz factor (gamma) ---
    # We can rearrange the energy conservation equation to solve for gamma.
    try:
        gamma_calculated = (2 * m_p_c2) / (4 * m_A_c2)
    except ZeroDivisionError:
        return "Constraint check failed: The mass of particle A (m_A) cannot be zero."

    # A physical constraint is that gamma must be >= 1.
    if gamma_calculated < 1:
        return (f"Constraint check failed: The calculated Lorentz factor gamma ({gamma_calculated:.4f}) is less than 1. "
                "This is physically impossible and implies the total rest mass of the final products "
                "is greater than the initial energy, which violates energy conservation.")

    # --- Step 3: Calculate the velocity (v) from gamma ---
    # The Lorentz factor is related to velocity by the formula: gamma = 1 / sqrt(1 - (v/c)^2).
    # We can solve for the ratio v/c: v/c = sqrt(1 - 1/gamma^2).
    v_over_c_calculated = np.sqrt(1 - 1/gamma_calculated**2)

    # --- Step 4: Compare the calculated result with the provided answer ---
    # We check if the calculated velocity is close to the velocity from option B (0.77c).
    # A tolerance is used to account for potential rounding in the options.
    tolerance = 0.005 

    if abs(v_over_c_calculated - provided_answer_v_over_c) < tolerance:
        return "Correct"
    else:
        return (f"Incorrect. The provided answer is v = {provided_answer_v_over_c}c (Option B), "
                f"but the correct calculation yields v/c = {v_over_c_calculated:.4f}. "
                f"The calculation steps are: E_initial = 2 * {m_p_c2:.2f} = {E_initial:.2f} MeV. "
                f"E_final = 4 * gamma * {m_A_c2:.2f} MeV. "
                f"Equating them gives gamma = {gamma_calculated:.4f}, which corresponds to v/c = {v_over_c_calculated:.4f}. "
                f"The provided answer of 0.77c is consistent with this result, but the check failed the tolerance. "
                f"Let's re-evaluate. The calculated value {v_over_c_calculated:.4f} rounds to 0.77. The answer is correct.")

# The logic in the final `else` block is a bit complex. Let's simplify.
# The calculation `0.7688...` rounds to `0.77`. So the answer is correct.
# The code should reflect this.

def final_check():
    m_p_c2 = 938.272
    m_A_c2 = 300.0
    provided_answer_v_over_c = 0.77
    
    # Calculate gamma
    gamma = (2 * m_p_c2) / (4 * m_A_c2)
    
    # Check physical constraint
    if gamma < 1:
        return f"Incorrect. Calculated gamma is {gamma}, which is less than 1 and physically impossible."
        
    # Calculate v/c
    v_over_c = np.sqrt(1 - 1/gamma**2)
    
    # Check if the provided answer is the closest option
    options = {'A': 0.86, 'B': 0.77, 'C': 0.91, 'D': 0.96}
    
    # Find the option that is numerically closest to our calculation
    closest_option_key = min(options, key=lambda k: abs(options[k] - v_over_c))
    
    if closest_option_key == 'B':
        # The answer B is indeed the closest and correct option.
        return "Correct"
    else:
        return (f"Incorrect. The provided answer is B (0.77c), but the calculation "
                f"yields v/c = {v_over_c:.4f}. This value is closest to option "
                f"{closest_option_key} ({options[closest_option_key]}c).")

# Execute the check
print(final_check())