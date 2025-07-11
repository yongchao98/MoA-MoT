import math

def solve_radiochemistry_problem():
    """
    This script solves the radiochemistry problem by calculating the time
    between chemical separation and the first measurement.
    """
    # Step 1: Define constants and identify the radioactive species.
    # The most likely pair is Ba-140/La-140 due to short irradiation and chemistry.
    T_Ba = 12.75  # Half-life of Ba-140 in days
    T_La = 40.28 / 24.0  # Half-life of La-140 in days (converted from hours)
    activity_ratio_measured = 2.1 / 1.4

    # Calculate decay constants (lambda = ln(2) / T_half) in units of day^-1
    lambda_Ba = math.log(2) / T_Ba
    lambda_La = math.log(2) / T_La

    print("--- Problem Analysis ---")
    print(f"Parent nuclide: Ba-140 (Half-life = {T_Ba:.2f} days, lambda = {lambda_Ba:.5f} day^-1)")
    print(f"Daughter nuclide: La-140 (Half-life = {T_La:.3f} days, lambda = {lambda_La:.5f} day^-1)")
    print("The measured Cherenkov activity is assumed to be from the daughter, La-140.")
    print("-" * 26)

    # Step 2: Define the function representing the ratio of activities.
    # The activity of La-140 at time 't' after separation is proportional to:
    # exp(-lambda_Ba * t) - exp(-lambda_La * t)
    def calculate_activity_ratio(t):
        """Calculates the ratio of La-140 activity at t+14 days to activity at t."""
        numerator = math.exp(-lambda_Ba * (t + 14)) - math.exp(-lambda_La * (t + 14))
        denominator = math.exp(-lambda_Ba * t) - math.exp(-lambda_La * t)
        if denominator == 0:
            return float('inf')
        return numerator / denominator

    # Step 3: Numerically solve for t1.
    # We search for the value of t1 that makes the calculated ratio equal to the measured ratio (1.5).
    t1_solution = 0
    min_difference = float('inf')

    # Perform a numerical search over a plausible time range (e.g., 0 to 5 days)
    for i in range(1, 5001):
        t_test = i * 0.001
        current_ratio = calculate_activity_ratio(t_test)
        difference = abs(current_ratio - activity_ratio_measured)
        if difference < min_difference:
            min_difference = difference
            t1_solution = t_test
    
    # Step 4: Print the equation and the final answer.
    print("\n--- Calculation ---")
    print("Let t1 be the time from separation to the first measurement.")
    print("We solve the following equation for t1:")
    print(f"{activity_ratio_measured:.1f} = [exp(-{lambda_Ba:.5f}*(t1+14)) - exp(-{lambda_La:.5f}*(t1+14))] / [exp(-{lambda_Ba:.5f}*t1) - exp(-{lambda_La:.5f}*t1)]")
    
    print("\n--- Result ---")
    print(f"The time from chemical separation to the first analysis (t1) is {t1_solution:.2f} days.")
    print("This is the 'age' of the prepared sample that can be determined from the data.")

solve_radiochemistry_problem()
<<<1.04>>>