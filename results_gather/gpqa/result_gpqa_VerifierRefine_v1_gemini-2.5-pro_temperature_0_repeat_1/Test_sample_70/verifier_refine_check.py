import math

def check_exoplanet_temperature_ratio():
    """
    This function verifies the answer to the exoplanet temperature ratio problem.

    The problem asks for the ratio of equilibrium temperatures between Planet_4 and Planet_2 (T4 / T2).
    The key relationships are:
    1. Equilibrium temperature (T_eq) is proportional to the inverse square root of the semi-major axis (a): T_eq ∝ a^(-1/2).
    2. Kepler's Third Law for planets around the same star states that the square of the orbital period (P) is proportional to the cube of the semi-major axis: P^2 ∝ a^3, which means a ∝ P^(2/3).

    Combining these, we get the relationship between temperature and period:
    T_eq ∝ (P^(2/3))^(-1/2) = P^(-1/3).

    Therefore, the ratio of temperatures is:
    T4 / T2 = (P4 / P2)^(-1/3) = (P2 / P4)^(1/3).
    """

    # Given orbital period ratios P1:P2:P3:P4:P5 = 1:2:2.5:3.5:5
    period_ratios = {
        "Planet_1": 1.0,
        "Planet_2": 2.0,
        "Planet_3": 2.5,
        "Planet_4": 3.5,
        "Planet_5": 5.0,
    }

    # Get the relative periods for Planet_2 and Planet_4
    p2 = period_ratios["Planet_2"]
    p4 = period_ratios["Planet_4"]

    # The answer to check is D, which corresponds to ~0.83
    answer_value = 0.83
    
    # Constraint check: The problem specifies 5 planets with given ratios.
    if len(period_ratios) != 5:
        return "Constraint failed: The problem specifies a system of five exoplanets."

    # Calculate the ratio of the periods P2 / P4
    period_ratio_p2_over_p4 = p2 / p4

    # Calculate the expected temperature ratio T4 / T2 = (P2 / P4)^(1/3)
    try:
        calculated_temp_ratio = period_ratio_p2_over_p4 ** (1/3)
    except (ValueError, TypeError):
        return "Calculation error: Could not compute the temperature ratio. Check input values."

    # Check if the calculated ratio is close to the provided answer's value
    # A tolerance of 0.01 is reasonable for a two-decimal-place answer.
    tolerance = 0.01
    if abs(calculated_temp_ratio - answer_value) < tolerance:
        return "Correct"
    else:
        return (f"Incorrect. The provided answer is {answer_value}, but the calculated value is approximately {calculated_temp_ratio:.4f}.\n"
                f"Reasoning:\n"
                f"1. The ratio of orbital periods P2/P4 is {p2}/{p4} = {period_ratio_p2_over_p4:.4f}.\n"
                f"2. The temperature ratio T4/T2 is derived as (P2/P4)^(1/3).\n"
                f"3. The correct value is ({period_ratio_p2_over_p4:.4f})^(1/3) ≈ {calculated_temp_ratio:.4f}, which rounds to {round(calculated_temp_ratio, 2)}.\n"
                f"The answer D (~0.83) is consistent with this calculation, but the provided value in the check ({answer_value}) might be slightly different if it were another option.")

# Run the check
print(check_exoplanet_temperature_ratio())