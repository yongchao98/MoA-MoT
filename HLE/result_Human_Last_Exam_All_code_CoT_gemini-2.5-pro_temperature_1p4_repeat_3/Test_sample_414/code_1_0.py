import math

def solve_radiochemistry_problem():
    """
    This script calculates the time between reactor irradiation and the first analysis
    of a radioactive sample based on two Cherenkov counting measurements.
    """
    # --- 1. Define constants and given values ---
    # Half-life of Yttrium-90 in hours
    T_half_Y90_h = 64.1
    # First measured activity in kBq/mL
    A1_activity = 1.4
    # Second measured activity in kBq/mL
    A2_activity = 2.1
    # Time for Y-90 ingrowth between measurements in days
    t_ingrowth_days = 14
    
    print("### Analysis of the Radiochemistry Problem ###\n")
    print("Step 1: Define physical constants and initial conditions.")
    print(f"Half-life of Y-90 (T_1/2): {T_half_Y90_h} hours")
    print(f"Activity at first measurement (A1): {A1_activity} kBq/mL")
    print(f"Activity at second measurement (A2): {A2_activity} kBq/mL after 14 days")
    print("-" * 50)

    # --- 2. Calculate decay constant and convert time ---
    # Calculate the decay constant for Y-90 (lambda)
    lambda_Y90_per_hour = math.log(2) / T_half_Y90_h
    # Convert ingrowth time from days to hours
    t_ingrowth_hours = t_ingrowth_days * 24

    # --- 3. Calculate the parent Sr-90 activity (A_Sr) ---
    print("Step 2: Calculate the parent Sr-90 activity (A_Sr) from the second measurement.")
    print("The second measurement (A2) was taken after 14 days of Y-90 ingrowth from a sample initially free of Y-90.")
    print("The formula for daughter ingrowth is: A_Y(t) = A_Sr * (1 - exp(-lambda * t))")
    print("We can rearrange this to solve for A_Sr:\n")
    
    # Calculation
    ingrowth_factor = 1 - math.exp(-lambda_Y90_per_hour * t_ingrowth_hours)
    A_Sr = A2_activity / ingrowth_factor
    
    print("Equation for Sr-90 Activity (A_Sr):")
    print(f"A_Sr = A2 / (1 - exp(-lambda_Y90 * t))")
    # Using format to show the numbers that are part of the equation
    print(f"A_Sr = {A2_activity} / (1 - exp(-({math.log(2):.4f} / {T_half_Y90_h}) * {t_ingrowth_hours}))")
    print(f"Resulting Sr-90 Activity: A_Sr = {A_Sr:.3f} kBq/mL")
    print("-" * 50)

    # --- 4. Calculate the initial cooling time (T_cool) ---
    print("Step 3: Calculate the cooling time (T_cool) from the first measurement.")
    print("The first measurement (A1) is the activity of Y-90 that grew in between irradiation and the separation.")
    print("Using the same ingrowth formula, with the now-known A_Sr, we solve for the time T_cool:\n")

    # Calculation
    try:
        ratio = A1_activity / A_Sr
        if ratio >= 1:
             print("Error: A1 activity cannot be greater than or equal to equilibrium activity A_Sr.")
             return
        T_cool_hours = -math.log(1 - ratio) / lambda_Y90_per_hour
    except ValueError:
        print("Error in calculation, check input values.")
        return

    print("Equation for Cooling Time (T_cool):")
    print("T_cool = -ln(1 - A1 / A_Sr) / lambda_Y90")
    print(f"T_cool = -ln(1 - {A1_activity} / {A_Sr:.3f}) / ({math.log(2):.4f} / {T_half_Y90_h})")
    print(f"Resulting Cooling Time: T_cool = {T_cool_hours:.1f} hours")
    print("-" * 50)
    
    print("\nFinal Answer:")
    print("The approximate time between sample irradiation and the first analysis is the calculated cooling time.")
    print(f"Time = {T_cool_hours:.1f} hours")
    
    # Store the final numerical answer for the platform
    global final_answer
    final_answer = T_cool_hours

# Execute the function
solve_radiochemistry_problem()