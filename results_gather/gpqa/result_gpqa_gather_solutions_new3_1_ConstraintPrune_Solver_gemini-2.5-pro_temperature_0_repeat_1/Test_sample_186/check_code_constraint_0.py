import math

def check_answer():
    """
    Checks the correctness of the answer by verifying the detectability of each star.
    """

    # --- Define Constants and Constraints ---

    # Latitude of Paranal Observatory (degrees)
    PARANAL_LATITUDE = -24.6

    # Limiting apparent magnitude for S/N=10 in 1hr, based on the answer's analysis
    # of ESPRESSO documentation.
    LIMITING_MAGNITUDE = 17.0

    # --- Star Data ---
    # A list of dictionaries, each representing a star with its properties.
    stars = [
        {
            "name": "a) Canopus",
            "declination": -52.7,  # degrees
            "apparent_magnitude": -0.74,
            "absolute_magnitude": None,
            "distance_pc": None
        },
        {
            "name": "b) Polaris",
            "declination": 89.26,  # degrees
            "apparent_magnitude": 1.98,
            "absolute_magnitude": None,
            "distance_pc": None
        },
        {
            "name": "c) Star at 10 pc",
            "declination": 0.0,
            "apparent_magnitude": None,
            "absolute_magnitude": 15.0,
            "distance_pc": 10.0
        },
        {
            "name": "d) Star at 200 pc",
            "declination": 0.0,
            "apparent_magnitude": None,
            "absolute_magnitude": 15.0,
            "distance_pc": 200.0
        },
        {
            "name": "e) Star at 5 pc",
            "declination": 0.0,
            "apparent_magnitude": None,
            "absolute_magnitude": 15.0,
            "distance_pc": 5.0
        },
        {
            "name": "f) Star at 50 pc",
            "declination": 0.0,
            "apparent_magnitude": None,
            "absolute_magnitude": 15.0,
            "distance_pc": 50.0
        }
    ]

    # --- Helper Functions ---

    def is_visible(declination, latitude):
        """Checks if a star is ever above the horizon at a given latitude."""
        # A star is visible if its maximum altitude is > 0.
        # Max altitude = 90 - |latitude - declination|
        max_altitude = 90 - abs(latitude - declination)
        return max_altitude > 0

    def calculate_apparent_magnitude(absolute_magnitude, distance_pc):
        """Calculates apparent magnitude using the distance modulus formula."""
        # m = M - 5 + 5 * log10(d)
        if distance_pc <= 0:
            return float('inf') # Not physically meaningful
        return absolute_magnitude - 5 + 5 * math.log10(distance_pc)

    # --- Main Verification Logic ---

    detectable_stars = []
    analysis_log = []

    for star in stars:
        # Step 1: Check visibility
        visible = is_visible(star["declination"], PARANAL_LATITUDE)
        if not visible:
            analysis_log.append(f"{star['name']}: Not Detectable. Reason: Not visible from Paranal (DEC={star['declination']}).")
            continue

        # Step 2: Determine apparent magnitude
        m_v = star["apparent_magnitude"]
        if m_v is None:
            m_v = calculate_apparent_magnitude(star["absolute_magnitude"], star["distance_pc"])

        # Step 3: Check brightness
        bright_enough = m_v <= LIMITING_MAGNITUDE
        if bright_enough:
            detectable_stars.append(star['name'])
            analysis_log.append(f"{star['name']}: Detectable. (Visible, m_V={m_v:.2f} <= {LIMITING_MAGNITUDE})")
        else:
            analysis_log.append(f"{star['name']}: Not Detectable. Reason: Too faint (m_V={m_v:.2f} > {LIMITING_MAGNITUDE}).")

    # --- Final Check ---
    
    # The provided answer concludes there are 3 detectable stars, which corresponds to option D.
    expected_count = 3
    expected_option = 'D'
    
    # Let's check if our code's conclusion matches the provided answer.
    calculated_count = len(detectable_stars)

    if calculated_count != expected_count:
        return (f"Incorrect. The code calculated {calculated_count} detectable stars, "
                f"but the answer states there are {expected_count}. "
                f"The detectable stars found were: {detectable_stars}.")

    # The mapping from count to option is: A=4, B=2, C=5, D=3
    option_map = {4: 'A', 2: 'B', 5: 'C', 3: 'D'}
    calculated_option = option_map.get(calculated_count)

    if calculated_option != expected_option:
        return (f"Incorrect. The calculated count of {calculated_count} corresponds to option {calculated_option}, "
                f"but the provided answer is {expected_option}.")

    # If all checks pass, the answer is consistent and correct.
    return "Correct"

# Run the check and print the result
result = check_answer()
print(result)