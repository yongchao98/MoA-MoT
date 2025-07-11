import math

def solve_cooling_time():
    """
    Calculates the time between sample irradiation and analysis based on
    the activities of Ba-140 and La-140.
    """
    # Step 1: Define known values from the problem.
    # Half-life of the parent, Ba-140 (in days)
    T_half_parent = 12.75
    # Half-life of the daughter, La-140 (in hours, then converted to days)
    T_half_daughter_hr = 40.28
    T_half_daughter = T_half_daughter_hr / 24.0

    # The problem implies an activity ratio can be determined from the measurements.
    # We hypothesize that the values 1.4 and 2.1 kBq/mL represent the activities
    # of the daughter (La-140) and parent (Ba-140) respectively, just before chemical separation.
    A_daughter = 1.4  # Activity of La-140 in kBq/mL
    A_parent = 2.1   # Activity of Ba-140 in kBq/mL
    
    # Step 2: Calculate the decay constants (lambda = ln(2) / T_half)
    lambda_p = math.log(2) / T_half_parent
    lambda_d = math.log(2) / T_half_daughter

    # Step 3: Calculate the activity ratio at the time of analysis
    activity_ratio = A_daughter / A_parent
    
    # Step 4: Use the equation for the activity ratio to solve for time 't'.
    # The equation is: R = [lambda_d / (lambda_d - lambda_p)] * [1 - exp(-(lambda_d - lambda_p) * t)]
    # We need to solve for 't'.

    # Rearranging the formula for t:
    # R / [lambda_d / (lambda_d - lambda_p)] = 1 - exp(-(lambda_d - lambda_p) * t)
    # exp(-(lambda_d - lambda_p) * t) = 1 - (R * (lambda_d - lambda_p) / lambda_d)
    # -(lambda_d - lambda_p) * t = ln(1 - (R * (lambda_d - lambda_p) / lambda_d))
    # t = -ln(1 - (R * (lambda_d - lambda_p) / lambda_d)) / (lambda_d - lambda_p)

    lambda_diff = lambda_d - lambda_p
    term_in_log = 1 - (activity_ratio * lambda_diff / lambda_d)
    
    # Check if the term inside the logarithm is positive
    if term_in_log <= 0:
        print("Error: The activity ratio is not physically possible for this decay chain.")
        return

    cooling_time_t = -math.log(term_in_log) / lambda_diff

    # Step 5: Print the results and the equation used.
    print("--- Radionuclide Properties ---")
    print(f"Parent (Ba-140) Half-Life: {T_half_parent:.2f} days")
    print(f"Daughter (La-140) Half-Life: {T_half_daughter:.2f} days")
    print(f"Parent Decay Constant (λ_p): {lambda_p:.4f} day^-1")
    print(f"Daughter Decay Constant (λ_d): {lambda_d:.4f} day^-1")
    print("\n--- Calculation Steps ---")
    print("Assuming the activities at the time of separation were:")
    print(f"Activity of Parent (Ba-140), A_p = {A_parent} kBq/mL")
    print(f"Activity of Daughter (La-140), A_d = {A_daughter} kBq/mL")
    print(f"Calculated Activity Ratio, R = A_d / A_p = {activity_ratio:.4f}")
    
    print("\nUsing the formula for radionuclide age determination:")
    print("R = [λ_d / (λ_d - λ_p)] * [1 - exp(-(λ_d - λ_p) * t)]")
    
    # Show the numbers in the final equation being solved
    lambda_ratio_term = lambda_d / lambda_diff
    print(f"{activity_ratio:.4f} = [{lambda_d:.4f} / ({lambda_d:.4f} - {lambda_p:.4f})] * [1 - exp(-({lambda_d:.4f} - {lambda_p:.4f}) * t)]")
    print(f"{activity_ratio:.4f} = {lambda_ratio_term:.4f} * [1 - exp(-{lambda_diff:.4f} * t)]")

    print("\n--- Final Answer ---")
    print(f"The approximate time between irradiation and the first analysis (cooling time) is: {cooling_time_t:.2f} days.")

solve_cooling_time()