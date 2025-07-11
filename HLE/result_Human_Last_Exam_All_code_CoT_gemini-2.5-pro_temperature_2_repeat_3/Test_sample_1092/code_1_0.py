import math

def calculate_fade_time():
    """
    Calculates the time in years for a Just Noticeable Fade (JNF) to occur
    for a light-sensitive object under specific conditions.
    """
    # --- Step 1: Define constants based on conservation science principles ---

    # The light dose required to cause one Just Noticeable Fade (JNF) for an
    # ISO Bluewool 1 rated material is approximately 50,000 lux-hours.
    jnf_dose = 50000  # units: lux-hours

    # --- Step 2: Define variables from the user's request and common assumptions ---

    # The light level the object is exposed to.
    lux_level = 50  # units: lux

    # A damage factor for "UV-rich" light. We assume it's twice as damaging
    # as a standard, low-UV light source.
    uv_damage_factor = 2

    # The number of hours the object is exposed per year.
    # A standard assumption for museum displays is ~3,000 hours/year.
    annual_exposure_hours = 3000  # units: hours/year

    # --- Step 3: Perform the calculation ---

    # Calculate the total effective light dose the object receives per year.
    annual_effective_dose = lux_level * uv_damage_factor * annual_exposure_hours

    # Calculate the number of years it will take to reach the JNF threshold.
    years_to_fade = jnf_dose / annual_effective_dose

    # --- Step 4: Print the results ---

    print("To find the time until the next Just Noticeable Fade, we use the following equation:")
    print("Years = (Total JNF Dose) / (Lux Level * UV Damage Factor * Annual Exposure Hours)\n")

    print("Using the provided values and standard assumptions:")
    # We use math.trunc() to show the integer values from the inputs/assumptions
    # for clarity in the equation being printed.
    print(f"Years = {math.trunc(jnf_dose)} / ({math.trunc(lux_level)} * {math.trunc(uv_damage_factor)} * {math.trunc(annual_exposure_hours)})")
    print(f"Years = {math.trunc(jnf_dose)} / {math.trunc(annual_effective_dose)}")

    print(f"\nResult: It will take {years_to_fade:.2f} years for the next just noticeable fade to occur.")

calculate_fade_time()