import math

def solve_titan_mass_problem():
    """
    Calculates the mass of a rock using Titan's 5-bit fractional arithmetic,
    finds the result with the minimum absolute error, and prints the derivation.
    """
    # Step 1: Define initial values as 5-bit fractions
    f_4_3 = (4, 3)
    f_pi = (22, 7)  # Best 5-bit approximation for pi
    f_r = (1, 2)    # r = 0.5 cm
    f_r3 = (f_r[0]**3, f_r[1]**3) # r^3 = (1/8)
    f_rho = (9, 10) # rho = 0.9 kg/cm3

    # Step 2: Formulate the unsimplified exact fraction
    # M = (4 * 22 * 1^3 * 9) / (3 * 7 * 2^3 * 10) = 792 / 1680
    # Simplifying this gives 33/70, but 33 and 70 are > 31.
    target_val = 33 / 70

    # Step 3: Find the best 5-bit fraction approximation for the target value.
    # A 5-bit unsigned integer can be from 0 to 31.
    limit = 31
    best_fraction = (0, 1)
    min_approx_error = float('inf')

    for d in range(1, limit + 1):
        # We want n/d ≈ target_val, so n ≈ target_val * d.
        # We round to the nearest integer for the best numerator.
        n = round(target_val * d)
        
        # Numerator must also be within the 5-bit limit.
        if n > limit:
            continue
        
        current_approx_error = abs(target_val - n / d)
        if current_approx_error < min_approx_error:
            min_approx_error = current_approx_error
            best_fraction = (int(n), int(d))

    # Step 4: The final calculated mass is the best 5-bit fraction
    m_calc_num, m_calc_den = best_fraction
    m_calc_val = m_calc_num / m_calc_den

    # Step 5: Calculate the true mass using high-precision values for error analysis
    m_true = (4/3) * math.pi * (0.5**3) * 0.9

    # Step 6: Calculate the final absolute error
    final_abs_error = abs(m_true - m_calc_val)
    
    # Step 7: Print the detailed derivation
    print("Derivation of the mass calculation:")
    print("Mass = (4/3) * pi * r^3 * rho")
    print(f"Using 5-bit fractional approximations: M = ({f_4_3[0]}/{f_4_3[1]}) * ({f_pi[0]}/{f_pi[1]}) * ({f_r3[0]}/{f_r3[1]}) * ({f_rho[0]}/{f_rho[1]})")
    print("Multiplying and simplifying these terms results in the fraction 33/70.")
    print("This is unrepresentable, as 33 and 70 are larger than 31.")
    print(f"We approximate 33/70 by finding the closest representable fraction n/d where n,d <= 31.")
    print(f"The best fractional approximation is {m_calc_num}/{m_calc_den}.")
    
    print("\nFinal Equation:")
    print(f"{m_calc_num} / {m_calc_den}")
    
    print(f"\nThis results in a calculated mass of {m_calc_val:.5f} kg.")
    print(f"The true mass is approximately {m_true:.5f} kg.")
    print(f"The absolute error is |{m_true:.5f} - {m_calc_val:.5f}| = {final_abs_error:.5f}")

    # The required answer 'e' is the smallest absolute error, rounded to 0.001.
    final_answer = round(final_abs_error, 3)
    print(f"\nThe smallest absolute error, rounded to 0.001, is: {final_answer}")
    print(f"<<<{final_answer}>>>")

solve_titan_mass_problem()