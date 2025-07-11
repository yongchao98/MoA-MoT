import math

def calculate_snr():
    """
    This script calculates the signal-to-noise ratio (SNR) at a specific location
    based on provided sound level measurements.
    """

    # --- Helper Functions ---
    def get_level_at_distance(source_level, distance):
        """Calculates sound level at a distance from a source."""
        if distance <= 0:
            return float('inf')
        return source_level - 20 * math.log10(distance)

    def get_source_level(level_at_distance, distance):
        """Calculates the source sound level (at 1m) from a measurement."""
        return level_at_distance + 20 * math.log10(distance)

    def subtract_decibels(total_db, component_db):
        """
        Subtracts a sound source (in dB) from a total combined sound level (in dB).
        Returns the level of the remaining sound source(s).
        """
        if component_db >= total_db:
            # If the component is louder than the total, the remainder is effectively zero (-inf dB)
            # This can also happen due to floating point inaccuracies if they are very close.
            # In our case, the component is much smaller.
            return -float('inf')

        total_intensity = 10**(total_db / 10)
        component_intensity = 10**(component_db / 10)
        remaining_intensity = total_intensity - component_intensity
        return 10 * math.log10(remaining_intensity)

    def combine_decibels(*db_levels):
        """Combines multiple sound levels into a single equivalent dB level."""
        total_intensity = sum(10**(db / 10) for db in db_levels)
        return 10 * math.log10(total_intensity)

    # --- Given Information ---

    # Signal is the sound of people talking at their location.
    signal_level = 75  # dB

    # Assumption: "talking at 75 dB" means the source level at 1m is 75 dB.
    people_source_level = 75  # dB at 1m

    # Measurement 1: Train and People
    meas1_total_level = 100  # dB
    meas1_dist_to_train = 10  # meters
    meas1_dist_to_people = 20  # meters

    # Measurement 2: Construction and People
    meas2_total_level = 115  # dB
    meas2_dist_to_construction = 20  # meters
    meas2_dist_to_people = 30  # meters

    # Distances to the final location (where the people are)
    dist_people_to_train = 30  # meters
    dist_people_to_construction = 50  # meters

    # --- Calculations ---

    # 1. Determine Train's Source Level
    people_level_at_meas1 = get_level_at_distance(people_source_level, meas1_dist_to_people)
    train_level_at_meas1 = subtract_decibels(meas1_total_level, people_level_at_meas1)
    train_source_level = get_source_level(train_level_at_meas1, meas1_dist_to_train)

    # 2. Determine Construction's Source Level
    people_level_at_meas2 = get_level_at_distance(people_source_level, meas2_dist_to_people)
    construction_level_at_meas2 = subtract_decibels(meas2_total_level, people_level_at_meas2)
    construction_source_level = get_source_level(construction_level_at_meas2, meas2_dist_to_construction)

    # 3. Calculate Noise Levels at the People's Location
    train_noise_at_location = get_level_at_distance(train_source_level, dist_people_to_train)
    construction_noise_at_location = get_level_at_distance(construction_source_level, dist_people_to_construction)

    # 4. Combine Noise Levels to get Total Noise
    total_noise_level = combine_decibels(train_noise_at_location, construction_noise_at_location)

    # 5. Calculate Final SNR
    snr = signal_level - total_noise_level

    # --- Output Final Result ---
    print("Final Equation: SNR = Signal Level - Total Noise Level")
    print(f"Signal Level (S): {signal_level:.2f} dB")
    print(f"Total Noise Level (N): {total_noise_level:.2f} dB")
    print(f"SNR = {signal_level:.2f} - {total_noise_level:.2f} = {snr:.2f} dB")

calculate_snr()