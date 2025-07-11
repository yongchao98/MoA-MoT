def calculate_fade_time():
    """
    Calculates the time in years for a just noticeable fade (JNF) to occur
    for a highly light-sensitive material under specific light conditions.
    """
    # Step 1: Standard exposure limit for a JNF on a highly sensitive material (ISO Bluewool 1).
    bw1_jnf_lux_hours = 50000.0

    # Step 2: Adjust for UV-rich light, which accelerates damage. We'll use a conservative factor of 2.
    uv_damage_factor = 2.0
    adjusted_jnf_lux_hours = bw1_jnf_lux_hours / uv_damage_factor

    # Step 3: Calculate the total annual light exposure.
    # Given light level in lux.
    lux_level = 50.0
    # Assuming a standard 8-hour "daily" exposure for a displayed object.
    hours_per_day = 8.0
    days_per_year = 365.0
    
    annual_exposure = lux_level * hours_per_day * days_per_year

    # Step 4: Calculate the time in years to reach the JNF threshold.
    years_to_fade = adjusted_jnf_lux_hours / annual_exposure

    print("The time to the next just noticeable fade is the total allowed exposure (adjusted for UV) divided by the annual exposure.")
    print("The calculation is:")
    print(f"({bw1_jnf_lux_hours} / {uv_damage_factor}) / ({lux_level} * {hours_per_day} * {days_per_year})")
    print(f"= {adjusted_jnf_lux_hours} / {annual_exposure}")
    print(f"= {years_to_fade:.2f} years")

calculate_fade_time()
<<<0.17>>>