import math

def calculate_total_spl():
    """
    Calculates the total sound pressure level (SPL) at a new location
    from multiple sound sources.
    """
    # 1. Define sound sources: name, SPL at 1m (dB), initial coordinates (x, y)
    sources = {
        'Dog': {'spl_at_1m': 55, 'coords': (-25, 0)},
        'Train': {'spl_at_1m': 110, 'coords': (50, 0)},
        'Construction': {'spl_at_1m': 90, 'coords': (0, 75)},
        'People': {'spl_at_1m': 75, 'coords': (0, -10)}
    }

    # 2. Define the listener's new position
    listener_new_pos = (0, 25)

    attenuated_spls = {}
    relative_intensities = {}

    print("--- Calculating SPL from each source at the new location ---\n")

    # Loop through each source to calculate new distance and attenuated SPL
    for name, data in sources.items():
        # 3. Calculate new distance
        source_coords = data['coords']
        distance = math.sqrt((source_coords[0] - listener_new_pos[0])**2 + (source_coords[1] - listener_new_pos[1])**2)

        # 4. Calculate attenuated SPL at the new distance
        spl_at_1m = data['spl_at_1m']
        # L2 = L1 - 20 * log10(r2)
        attenuated_spl = spl_at_1m - 20 * math.log10(distance)
        attenuated_spls[name] = attenuated_spl

        # 5. Convert attenuated SPL to relative intensity
        # I = 10^(L/10)
        relative_intensities[name] = 10**(attenuated_spl / 10)

        print(f"Source: {name}")
        print(f"  New distance: {distance:.2f} meters")
        print(f"  Attenuated SPL: {attenuated_spl:.2f} dB\n")

    # 6. Combine the sound levels
    # Sum the relative intensities
    total_relative_intensity = sum(relative_intensities.values())

    # Convert the total intensity back to decibels
    total_spl = 10 * math.log10(total_relative_intensity)

    # --- Final Output ---
    print("--- Combining Sound Levels ---")
    print("To find the total sound level, we convert each dB value to a linear relative intensity,")
    print("sum them, and then convert back to decibels.")
    print("Formula: Total dB = 10 * log10( sum(10^(L_i / 10)) )\n")

    print("The final calculation is:")
    
    # Building the string for the equation
    intensity_values = list(relative_intensities.values())
    sum_string = " + ".join([f"10^({attenuated_spls[name]:.2f}/10)" for name in sources.keys()])
    numeric_sum_string = " + ".join([f"{val:.2f}" for val in intensity_values])
    
    print(f"Total dB = 10 * log10( {sum_string} )")
    print(f"Total dB = 10 * log10( {numeric_sum_string} )")
    print(f"Total dB = 10 * log10( {total_relative_intensity:.2f} )")
    print(f"\nFinal Total Sound Level: {total_spl:.2f} dB")


if __name__ == "__main__":
    calculate_total_spl()
    # The calculated value is ~75.1 dB. We'll use a precise calculation for the final answer.
    # Recalculating here just for the final tag, using values from execution.
    # d_dog = 35.355, spl_dog = 24.03
    # d_train = 55.901, spl_train = 75.05
    # d_const = 50.0, spl_const = 56.02
    # d_people = 35.0, spl_people = 44.12
    # i_dog = 10**(2.403) = 252.9
    # i_train = 10**(7.505) = 31,989,016.5
    # i_const = 10**(5.602) = 399,945.1
    # i_people = 10**(4.412) = 25,822.5
    # sum = 32,415,037
    # 10*log10(sum) = 75.107
    final_answer = 75.11 # Rounded to two decimal places
    # print(f"\n<<<{final_answer}>>>")