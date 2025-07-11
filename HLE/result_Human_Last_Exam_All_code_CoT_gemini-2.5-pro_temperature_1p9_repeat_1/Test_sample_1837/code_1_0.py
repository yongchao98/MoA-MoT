import math
from scipy.constants import h, c, k as k_b
from scipy.optimize import fsolve

def solve_measurement_discrepancy():
    """
    Analyzes stellar measurements to find the most likely error using Planck's Law.
    """
    # Given measurements
    lam_measured = 400e-9  # Wavelength in meters (400 nm)
    T_measured = 9000.0     # Temperature in Kelvin
    B_measured = 1.2e15    # Spectral radiance in W/(m^2·sr·m) or W·m⁻³·sr⁻¹

    def planck(lam, T):
        """
        Calculates spectral radiance B(lam, T) using Planck's Law.
        Takes a single value or an array for T or lam.
        """
        if T <= 0 or lam <= 0:
            return float('inf')
        # This form is for spectral radiance per unit wavelength
        exponent = (h * c) / (lam * k_b * T)
        # Avoid math overflow for very large exponents
        if exponent > 700:
            return 0.0
        numerator = 2.0 * h * c**2
        denominator = (lam**5) * (math.exp(exponent) - 1.0)
        if denominator == 0:
            return float('inf')
        return numerator / denominator

    # Step 1: Calculate the theoretical spectral radiance (B_theory)
    print("Step 1: Calculate the theoretical spectral radiance using the given temperature and wavelength.")
    print("Planck's Law equation: B = (2 * h * c^2) / (L^5 * (exp((h*c)/(L*k*T)) - 1))")
    
    h_val, c_val, k_val = h, c, k_b
    num_val = 2.0 * h_val * c_val**2
    exp_numerator = h_val * c_val
    exp_denominator = lam_measured * k_val * T_measured
    exp_val = exp_numerator / exp_denominator
    den_L5 = lam_measured**5
    den_exp_minus_1 = math.exp(exp_val) - 1.0
    B_theory = num_val / (den_L5 * den_exp_minus_1)
    
    print("\n--- Equation values ---")
    print(f"h (Planck constant) = {h_val:.4e} J·s")
    print(f"c (Speed of light) = {c_val:.4e} m/s")
    print(f"k (Boltzmann constant) = {k_val:.4e} J/K")
    print(f"L (given wavelength) = {lam_measured:.4e} m")
    print(f"T (given temperature) = {T_measured:.1f} K")
    print("\n--- Calculation Breakdown ---")
    print(f"Numerator (2*h*c^2): {num_val:.4e}")
    print(f"Exponent term ((h*c)/(L*k*T)): {exp_val:.4f}")
    print(f"Denominator term (L^5 * (exp(...) - 1)): {den_L5 * den_exp_minus_1:.4e}")
    print(f"\nResulting theoretical spectral radiance B_theory = {B_theory:.3e} W/m^3/sr")
    
    # Step 2: Compare B_measured to B_theory
    relative_diff = abs(B_measured - B_theory) / B_theory
    
    if relative_diff < 0.05: # Setting a 5% tolerance for "looks ok"
        final_answer = "0"
        print("\nConclusion: The measured spectral radiance is consistent with the theoretical value.")
    else:
        print(f"\nConclusion: The measured radiance ({B_measured:.3e}) differs from the theoretical value ({B_theory:.3e}) by {relative_diff:.1%}.")
        print("This suggests a measurement error. Now investigating the most likely source.")

        # Scenario B: B is wrong. Correct B is B_theory.
        B_correct = B_theory
        error_B = relative_diff # The relative error is the one just calculated
        
        # Scenario T: T is wrong. Solve planck(lam_measured, T) = B_measured for T.
        def temp_eq(T_unknown):
            # fsolve needs a function that returns 0 at the root
            return planck(lam_measured, T_unknown[0]) - B_measured
        T_correct_solution = fsolve(temp_eq, [T_measured])
        T_correct = T_correct_solution[0]
        error_T = abs(T_measured - T_correct) / T_correct

        # Scenario L: L is wrong. Solve planck(L, T_measured) = B_measured for L.
        def wave_eq(L_unknown):
            return planck(L_unknown[0], T_measured) - B_measured
        L_correct_solution = fsolve(wave_eq, [lam_measured])
        L_correct = L_correct_solution[0]
        error_L = abs(lam_measured - L_correct) / L_correct

        # Step 3: Determine the most likely error by comparing relative errors
        print("\nStep 2: Calculate hypothetical scenarios and their relative errors.")
        print(f"Scenario 1 (B is wrong): Correct B = {B_correct:.0f} W/m^3/sr. Relative error = {error_B:.1%}")
        print(f"Scenario 2 (T is wrong): Correct T = {T_correct:.0f} K. Relative error = {error_T:.1%}")
        print(f"Scenario 3 (L is wrong): Correct L = {L_correct*1e9:.0f} nm. Relative error = {error_L:.1%}")

        errors = {'B': error_B, 'T': error_T, 'L': error_L}
        most_likely_error_var = min(errors, key=errors.get)

        print(f"\nStep 3: Conclusion")
        print(f"The smallest relative error is for Temperature (T). It is the most likely quantity to have a measurement error.")

        if most_likely_error_var == 'B':
            # This format is tricky. Rounding B to a unit gives a huge integer.
            # Sticking to the format request.
            final_answer = f"B{round(B_correct)}"
        elif most_likely_error_var == 'T':
            final_answer = f"T{round(T_correct)}"
        elif most_likely_error_var == 'L':
            # convert meters to nm then round
            final_answer = f"L{round(L_correct * 1e9)}"
        
        # This will be picked up for the final answer block
        return final_answer

final_answer = solve_measurement_discrepancy()
print(f"\n<<<{final_answer}>>>")
