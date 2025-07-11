def find_most_negative_factor():
    """
    Analyzes meteorological data to find the value most detrimental to tornadogenesis.

    The provided data shows an environment with powerful ingredients for severe storms:
    - High Instability: MLCAPE is 1567 J/kg and MUCAPE is 2136 J/kg.
    - Ample Moisture: Precipitable water is 1.4 inches and the SFC LCL is very low at 228m.
    - Extreme Wind Shear: SFC-1km SRH is 361 m^2/s^2 and SFC-6km shear is 66 kts.
    - High Composite Values: The Significant Tornado Parameter (STP) is 4.7.

    Despite these overwhelmingly favorable conditions, there is one kinematic value that stands out as a significant negative factor.

    The 'Corfidi Downshear' vector is 242/96 kt. The speed component of this vector is 96 kts.
    This indicates an extremely fast storm motion. While a result of very strong winds aloft,
    a storm moving this quickly (approximately 110 mph or 178 km/h) can struggle to
    produce and maintain a tornado. The storm's updraft may move too rapidly for the
    low-level rotation to properly organize, tighten, and be maintained at the surface.
    This kinematic extremity is the most significant factor that would inhibit or reduce
    the likelihood of an efficient, long-lived tornado in this environment.
    """
    negative_factor_name = "Corfidi Downshear Speed (knots)"
    negative_factor_value = 96
    
    print(f"The parameter that has the greatest negative impact on tornadogenesis is the {negative_factor_name}.")
    print(f"Its value is: {negative_factor_value}")

find_most_negative_factor()
