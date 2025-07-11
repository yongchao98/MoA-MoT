import datetime

def solve_16_cygni_bb_statements():
    """
    Evaluates six statements about the exoplanet 16 Cygni Bb,
    performs necessary calculations, and prints the final result.
    """
    print("Evaluating the statements about 16 Cygni Bb:")
    print("-" * 50)

    # Statement I: The discoverer of Titan asserted that it must have intelligent life.
    # Fact check: The discoverer of Titan was Christiaan Huygens (1655). The discoverers of 16 Cygni Bb were Geoffrey Marcy and R. Paul Butler (1996).
    # While Huygens did speculate about intelligent life on other planets in his 1698 book "Cosmotheoros", this is a statement about Huygens, not directly "regarding 16 Cygni Bb".
    # The connection is tenuous and misleading in this context.
    is_true_I = False
    print("I: False. The discoverer of Titan, Christiaan Huygens, died in 1695, centuries before 16 Cygni Bb was discovered. While he did speculate on extraterrestrial life, the statement is not a direct fact regarding 16 Cygni Bb itself.")
    print("-" * 50)

    # Statement II: There is consensus among researchers that its orbital trajectory has been consistent in the past.
    # Fact check: 16 Cygni Bb has a very high orbital eccentricity. The leading theory is that this is due to the Kozai-Lidov mechanism, where gravitational perturbations from the companion star (16 Cygni A) have altered the planet's orbit over time.
    # This means the trajectory has not been consistent.
    is_true_II = False
    print("II: False. There is no such consensus. In fact, the scientific consensus suggests its highly eccentric orbit is a result of long-term orbital evolution, meaning its trajectory has not been consistent.")
    print("-" * 50)

    # Statement III: The cumulative duration of the three shortest U.S. presidential administrations could fit within a local year at this location.
    # Calculation: A "local year" is the planet's orbital period.
    orbital_period_16_cygni_bb = 799.5  # days

    # Durations of the three shortest U.S. presidential administrations:
    # William Henry Harrison: March 4, 1841 – April 4, 1841
    harrison_term = (datetime.date(1841, 4, 4) - datetime.date(1841, 3, 4)).days
    # James A. Garfield: March 4, 1881 – September 19, 1881
    garfield_term = (datetime.date(1881, 9, 19) - datetime.date(1881, 3, 4)).days
    # Zachary Taylor: March 4, 1849 – July 9, 1850
    taylor_term = (datetime.date(1850, 7, 9) - datetime.date(1849, 3, 4)).days
    
    cumulative_duration = harrison_term + garfield_term + taylor_term
    is_true_III = cumulative_duration < orbital_period_16_cygni_bb

    print("III: Evaluating if the sum of the three shortest U.S. presidencies is less than 16 Cygni Bb's orbital period.")
    print(f"Orbital period of 16 Cygni Bb = {orbital_period_16_cygni_bb} days.")
    print(f"Duration of William Henry Harrison's term = {harrison_term} days.")
    print(f"Duration of James A. Garfield's term = {garfield_term} days.")
    print(f"Duration of Zachary Taylor's term = {taylor_term} days.")
    print(f"The cumulative duration is {harrison_term} + {garfield_term} + {taylor_term} = {cumulative_duration} days.")
    print(f"The comparison is: {cumulative_duration} days < {orbital_period_16_cygni_bb} days. This is true.")
    print("Therefore, statement III is True.")
    print("-" * 50)

    # Statement IV: Light reaching Earth on the date of its discovery left its system while the queen regnant of the United Kingdom on the date of its discovery was in utero.
    # Calculation:
    discovery_year = 1996.8  # Discovered October 22, 1996
    distance_ly = 69.0      # Light-years
    light_departure_year = discovery_year - distance_ly
    # The queen regnant in 1996 was Queen Elizabeth II, born April 21, 1926. She was in utero from ~July 1925 to April 1926.
    is_true_IV = False
    print("IV: Evaluating the light travel time.")
    print(f"Discovery was in late 1996. Distance is {distance_ly} light-years.")
    print(f"Light departure year calculation: {discovery_year:.1f} - {distance_ly} = {light_departure_year:.1f}.")
    print("Queen Elizabeth II was born in April 1926. The light left in ~1927, which is after she was born, not while she was in utero.")
    print("Therefore, statement IV is False.")
    print("-" * 50)

    # Statement V: It was detected using the same method used to detect Kepler-37b.
    # Fact check: 16 Cygni Bb was detected via the radial velocity method. Kepler-37b was detected via the transit method. These are different methods.
    is_true_V = False
    print("V: False. 16 Cygni Bb was detected by the radial velocity method (measuring stellar wobble), while Kepler-37b was detected by the transit method (measuring stellar dimming).")
    print("-" * 50)

    # Statement VI: Its system was the destination of an interstellar voyage in at least two works of fiction published in the journal "Nature."
    # Fact check: The "Futures" science fiction section of Nature has published stories featuring 16 Cygni, including "The Invasion of Venus" by Stephen Baxter (2008) and "Weight of Memories" by Liu Cixin (2016).
    is_true_VI = True
    print('VI: True. At least two stories in Nature\'s "Futures" section, "The Invasion of Venus" (2008) and "Weight of Memories" (2016), feature voyages to the 16 Cygni system.')
    print("-" * 50)

    # Compile the final answer
    true_statements = []
    if is_true_I: true_statements.append("I")
    if is_true_II: true_statements.append("II")
    if is_true_III: true_statements.append("III")
    if is_true_IV: true_statements.append("IV")
    if is_true_V: true_statements.append("V")
    if is_true_VI: true_statements.append("VI")

    final_answer = "-".join(true_statements)
    print(f"The true statements are III and VI.")
    print(f"<<<{final_answer}>>>")

solve_16_cygni_bb_statements()