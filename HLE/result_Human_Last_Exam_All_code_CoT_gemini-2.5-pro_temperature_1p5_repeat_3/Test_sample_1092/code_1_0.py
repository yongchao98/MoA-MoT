import math

def calculate_fade_time():
    """
    Calculates the time in years for a Just Noticeable Fade (JNF) to occur
    on a light-sensitive object based on its ISO Bluewool Rating and exposure conditions.
    """

    # --- Constants and Assumptions ---

    # JNF threshold for ISO Bluewool 1 is ~50,000 lux-hours.
    JNF_EXPOSURE_BW1 = 50000.0

    # Given light level in lux.
    LUX_LEVEL = 50.0

    # Assumption: "UV-rich" light is twice as damaging as standard UV-filtered light.
    UV_DAMAGE_FACTOR = 2.0

    # Assumption: "Daily" exposure means 8 hours per day.
    HOURS_PER_DAY = 8.0

    # Average number of days in a year.
    DAYS_PER_YEAR = 365.25

    # --- Calculation ---

    # Calculate the total annual light exposure, adjusted for UV damage.
    # Annual Exposure = Lux * UV Factor * Hours/Day * Days/Year
    annual_exposure = LUX_LEVEL * UV_DAMAGE_FACTOR * HOURS_PER_DAY * DAYS_PER_YEAR

    # Calculate the number of years to reach the JNF threshold.
    # Years to Fade = Total JNF Exposure / Annual Exposure
    years_to_fade = JNF_EXPOSURE_BW1 / annual_exposure

    # --- Output ---

    print("To calculate the time to a Just Noticeable Fade (JNF), we use the following formula:")
    print("Years = (Total Exposure for JNF) / (Annual Exposure)\n")

    print("Here are the values used in the calculation:")
    print(f"- Total Exposure for JNF (Bluewool 1): {int(JNF_EXPOSURE_BW1)} lux-hours")
    print(f"- Light Level: {int(LUX_LEVEL)} lux")
    print(f"- Assumed UV Damage Factor: {int(UV_DAMAGE_FACTOR)}")
    print(f"- Assumed Hours of Exposure per Day: {int(HOURS_PER_DAY)}")
    print(f"- Days per Year: {DAYS_PER_YEAR}\n")

    print("The final calculation is:")
    print(f"Years to fade = {int(JNF_EXPOSURE_BW1)} / ({int(LUX_LEVEL)} * {int(UV_DAMAGE_FACTOR)} * {int(HOURS_PER_DAY)} * {DAYS_PER_YEAR})")
    print(f"Result: {years_to_fade:.2f} years")


if __name__ == "__main__":
    calculate_fade_time()
<<<0.17>>>