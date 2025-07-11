def calculate_fade_time():
    """
    Calculates the time in years for a just noticeable fade (JNF) to occur
    on a highly sensitive material under specific light conditions.
    """
    # Step 1 & 2: Define constants and the UV-adjusted exposure for a JNF.
    # ISO Bluewool 1 requires ~50,000 lux-hours for a JNF under UV-filtered light.
    base_exposure_for_jnf = 50000.0  # lux-hours
    # UV-rich light is estimated to be ~3 times more damaging.
    uv_damage_factor = 3.0
    # The light dose for one JNF under UV-rich light.
    adjusted_exposure_for_jnf = base_exposure_for_jnf / uv_damage_factor

    # Step 3: Calculate the annual light exposure.
    illuminance_lux = 50.0
    # Assuming a standard 8-hour exhibition day.
    hours_per_day = 8.0
    days_per_year = 365.0
    annual_exposure = illuminance_lux * hours_per_day * days_per_year

    # Step 4: Calculate the time in years to reach one JNF.
    years_to_jnf = adjusted_exposure_for_jnf / annual_exposure

    # Output the explanation, the full equation, and the result.
    print("The time to the next Just Noticeable Fade (JNF) is calculated as follows:")
    print("\n1. Determine the light dose for one JNF, adjusted for UV-rich light:")
    print(f"   Adjusted Dose = Base Exposure for JNF / UV Damage Factor")
    print(f"   Adjusted Dose = {base_exposure_for_jnf} / {uv_damage_factor} = {adjusted_exposure_for_jnf:.2f} lux-hours")

    print("\n2. Calculate the total annual light exposure:")
    print(f"   Annual Exposure = Illuminance * Hours per Day * Days per Year")
    print(f"   Annual Exposure = {illuminance_lux} * {hours_per_day} * {days_per_year} = {annual_exposure} lux-hours/year")

    print("\n3. Calculate the years until the next JNF:")
    print(f"   Years = Adjusted Dose / Annual Exposure")
    print(f"   Years = {adjusted_exposure_for_jnf:.2f} / {annual_exposure}")
    print(f"\nResult: The next just noticeable fade will occur in approximately {years_to_jnf:.2f} years.")

calculate_fade_time()