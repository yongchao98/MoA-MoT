import math

def calculate_fading_time():
    """
    Calculates the time in years for a just noticeable fade to occur on a
    light-sensitive object.
    """
    # Step 1: Define the parameters based on the problem and conservation standards.
    # The total light exposure limit for a Just Noticeable Fade (JNF)
    # for an ISO Bluewool 1 rated material.
    total_lux_hours_for_jnf = 50000

    # The intensity of the light exposure in lux.
    lux_level = 50

    # Assuming a standard exhibition time of 8 hours per day.
    hours_per_day = 8
    
    # Using 365.25 to average in leap years for better accuracy over time.
    days_per_year = 365.25

    # Step 2: Calculate the total annual light exposure.
    # Annual Exposure = Lux Level * Hours per Day * Days per Year
    annual_exposure = lux_level * hours_per_day * days_per_year

    # Step 3: Calculate the years until fading occurs.
    # Years to Fade = Total Fading Threshold / Annual Exposure
    years_to_fade = total_lux_hours_for_jnf / annual_exposure

    # Step 4: Print the equation and the final answer.
    print("The time to fading is calculated by dividing the material's total exposure limit by its annual exposure.")
    print("\nFinal Equation:")
    print(f"Years to Fade = {total_lux_hours_for_jnf} / ({lux_level} * {hours_per_day} * {days_per_year})")
    print(f"Result: The object will show a just noticeable fade in approximately {years_to_fade:.3f} years.")

# Run the calculation
calculate_fading_time()
