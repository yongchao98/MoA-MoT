import math

def solve_cygni_statements():
    """
    Analyzes six statements about 16 Cygni Bb, prints the reasoning for each,
    and determines the final sequence of true statements.
    """
    print("Analyzing the statements about 16 Cygni Bb:\n")
    
    true_statements = []

    # --- Statement I ---
    # The discoverer of Titan asserted that it must have intelligent life.
    print("--- Checking Statement I ---")
    print("Fact: The discoverer of Titan was Christiaan Huygens (lived 1629-1695).")
    print("Fact: 16 Cygni Bb was discovered in 1996.")
    print("Conclusion: Huygens died over 300 years before 16 Cygni Bb was discovered, so he could not have made any assertion about it.")
    print("Statement I is FALSE.\n")

    # --- Statement II ---
    # There is consensus among researchers that its orbital trajectory has been consistent in the past.
    print("--- Checking Statement II ---")
    print("Fact: 16 Cygni is a triple-star system. The gravity of the third star, 16 Cygni A, perturbs the planet's orbit around 16 Cygni B.")
    print("Fact: The prevailing theory is the Kozai-Lidov mechanism, which suggests the planet's orbital eccentricity and inclination vary significantly over long periods.")
    print("Conclusion: The consensus is that the orbit has NOT been consistent.")
    print("Statement II is FALSE.\n")

    # --- Statement III ---
    # The cumulative duration of the three shortest U.S. presidential administrations could fit within a local year at this location.
    print("--- Checking Statement III ---")
    local_year_days = 799.5
    print(f"A local year on 16 Cygni Bb (its orbital period) is {local_year_days} Earth days.")
    # Durations of the three shortest U.S. presidential administrations
    harrison_days = 31  # William Henry Harrison
    garfield_days = 199 # James A. Garfield
    taylor_days = 492   # Zachary Taylor (1 year, 127 days)
    print(f"The three shortest U.S. presidential administrations lasted:")
    print(f"- William Henry Harrison: {harrison_days} days")
    print(f"- James A. Garfield: {garfield_days} days")
    print(f"- Zachary Taylor: {taylor_days} days")
    total_duration_days = harrison_days + garfield_days + taylor_days
    print(f"Equation: {harrison_days} + {garfield_days} + {taylor_days} = {total_duration_days} days")
    is_shorter = total_duration_days < local_year_days
    print(f"Comparison: Is {total_duration_days} days < {local_year_days} days? {is_shorter}")
    if is_shorter:
        print("Conclusion: The cumulative duration fits within a local year.")
        print("Statement III is TRUE.\n")
        true_statements.append("III")
    else:
        print("Conclusion: The cumulative duration does not fit within a local year.")
        print("Statement III is FALSE.\n")

    # --- Statement IV ---
    # Light reaching Earth on the date of its discovery left its system while the queen regnant of the
    # United Kingdom on the date of its discovery was in utero.
    print("--- Checking Statement IV ---")
    discovery_year = 1996
    distance_ly = 69.0 # Distance in light-years
    queen_birth_year = 1926 # Queen Elizabeth II's birth year
    print(f"The planet was discovered in {discovery_year}.")
    print(f"The queen regnant at the time was Queen Elizabeth II, born in {queen_birth_year}.")
    print(f"The system is approximately {distance_ly} light-years away.")
    light_departure_year = discovery_year - distance_ly
    print(f"Equation: {discovery_year} - {distance_ly} = {int(light_departure_year)}")
    print(f"The light that reached Earth in {discovery_year} left its system in {int(light_departure_year)} (i.e., 1927).")
    # Queen Elizabeth II was in utero from approx. July 1925 to April 1926.
    print(f"Conclusion: In {int(light_departure_year)}, Queen Elizabeth II had already been born; she was not in utero.")
    print("Statement IV is FALSE.\n")

    # --- Statement V ---
    # It was detected using the same method used to detect Kepler-37b.
    print("--- Checking Statement V ---")
    print("Fact: 16 Cygni Bb was detected by the radial velocity method (observing stellar wobble).")
    print("Fact: Kepler-37b was detected by the transit method (observing stellar dimming).")
    print("Conclusion: The detection methods are different.")
    print("Statement V is FALSE.\n")
    
    # --- Statement VI ---
    # Its system was the destination of an interstellar voyage in at least two works of fiction published in the journal "Nature."
    print("--- Checking Statement VI ---")
    print("Fact: The journal Nature's 'Futures' section has published multiple sci-fi stories set in this system, including:")
    print("1. 'Worm's-Eye View' by Ken Brady (Nature, 2004)")
    print("2. 'A Letter from the Horse-Latitudes' by Robert J. Scherrer (Nature, 2005)")
    print("Conclusion: The system was a destination in at least two such works.")
    print("Statement VI is TRUE.\n")
    true_statements.append("VI")

    # --- Final Answer ---
    print("--- FINAL RESULT ---")
    final_answer = "-".join(true_statements)
    print(f"The true statements are III and VI.")
    print(f"Formatted as requested: {final_answer}")


if __name__ == "__main__":
    solve_cygni_statements()