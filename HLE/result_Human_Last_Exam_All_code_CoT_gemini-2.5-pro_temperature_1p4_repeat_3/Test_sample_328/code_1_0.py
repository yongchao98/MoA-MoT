import math

def estimate_mating_age():
    """
    This function estimates the average age of a mating male for two insect species
    based on their life history traits.
    """

    # --- Species 1: Tridactylophagus tartari (Scelionid Wasp) ---
    species_1_name = "Tridactylophagus tartari"

    # Life history traits:
    # These wasps mature very quickly after eclosion.
    s1_maturity_days = 1
    # They have a very short adult lifespan, typically just several days. We'll use 5 as an average.
    s1_lifespan_days = 5

    # The average age of a mating male is the midpoint of its reproductive period.
    # Reproductive period starts at maturity (day 1) and ends at death (day 5).
    s1_avg_age = (s1_maturity_days + s1_lifespan_days) / 2
    # The result is 3, which is a plausible average age in days.
    s1_estimate = "three days"


    # --- Species 2: Periplaneta americana (American Cockroach) ---
    species_2_name = "Periplaneta americana"

    # Life history traits:
    # Cockroaches mature relatively quickly, but have a long adult lifespan.
    # Maturity is about 1 week, which is negligible compared to the total lifespan for this estimate.
    s2_maturity_months = 0.25
    # Adult lifespan is approximately 1 year.
    s2_lifespan_months = 12

    # The average age of a mating male is the midpoint of its reproductive lifespan.
    # Since maturity is very short compared to the lifespan, we can approximate the average
    # age as half the total adult lifespan.
    # The average age would be roughly in the middle of their 1-year reproductive life.
    s2_avg_age = (s2_maturity_months + s2_lifespan_months) / 2
    # The result is ~6.125, which rounds nicely to 6 months.
    s2_estimate = "six months"


    # --- Print the results ---
    print("Estimating the average age of a mating male since eclosion (becoming an adult):")
    print("-" * 60)

    print(f"Species: {species_1_name}")
    print(f"This wasp has a short adult lifespan (avg. {s1_lifespan_days} days) and matures quickly (day {s1_maturity_days}).")
    print(f"Average age calculation: ({s1_maturity_days} + {s1_lifespan_days}) / 2 = {s1_avg_age:.0f} days.")
    print(f"The best average estimate is {s1_estimate}.")
    print("-" * 60)

    print(f"Species: {species_2_name}")
    print(f"This cockroach has a long adult lifespan (avg. {s2_lifespan_months} months).")
    print("An observed mating can happen anytime during this period.")
    print(f"Average age calculation (in months): (0 + {s2_lifespan_months}) / 2 = {math.ceil(s2_avg_age)} months.")
    print(f"The best average estimate is {s2_estimate}.")
    print("-" * 60)

    print(f"Final combined estimate: {s1_estimate} for {species_1_name}, {s2_estimate} for {species_2_name}.")

estimate_mating_age()
<<<C>>>