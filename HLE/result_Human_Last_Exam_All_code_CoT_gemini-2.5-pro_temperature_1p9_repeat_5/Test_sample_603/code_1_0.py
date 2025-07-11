# Plan:
# 1. The core limitation of bulk experiments is that they measure the average behavior of a massive number of molecules, obscuring any underlying differences, or "heterogeneity," within the sample.
# 2. I will simulate a sample containing two distinct populations of DNA molecules, each with its own melting temperature (Tm).
# 3. I will then calculate what a bulk measurement would observe by averaging the responses of these two populations.
# 4. The output will show that the bulk measurement provides a single, averaged value that does not reveal the presence of the two separate populations. This directly illustrates the concept of "impossibility to capture heterogeneity."

def simulate_fraction_denatured(temperature, melting_temp):
    """
    A simple sigmoid function to model the fraction of denatured molecules.
    A more realistic model would use thermodynamics (van 't Hoff equation),
    but this is sufficient to illustrate the concept.
    """
    steepness = 0.5  # Controls how sharp the transition is
    fraction = 1 / (1 + 2.71828**(-steepness * (temperature - melting_temp)))
    return fraction

# Let's define two distinct molecule populations with different melting points (Tm)
population_1_tm = 75.0  # e.g., perfectly formed duplexes
population_2_tm = 82.0  # e.g., duplexes with a stabilizing modification

# Let's pick a temperature to probe the system, for example, 78 degrees C
test_temperature = 78.0

# In a single-molecule experiment, we could isolate and measure each population
denatured_pop1 = simulate_fraction_denatured(test_temperature, population_1_tm)
denatured_pop2 = simulate_fraction_denatured(test_temperature, population_2_tm)

# In a bulk experiment, we measure the average of the whole sample.
# Assuming equal amounts of each population for simplicity.
bulk_measurement = (denatured_pop1 + denatured_pop2) / 2

print("This simulation shows how bulk measurements hide underlying molecular differences.\n")
print(f"Our sample has two distinct populations with melting temperatures (Tm) of {population_1_tm}°C and {population_2_tm}°C.")
print("-" * 60)
print(f"At a test temperature of {test_temperature}°C:")
print(f"   - Population 1 is {denatured_pop1:.2f} (or {denatured_pop1*100:.0f}%) denatured.")
print(f"   - Population 2 is {denatured_pop2:.2f} (or {denatured_pop2*100:.0f}%) denatured.")
print("-" * 60)
print("A bulk experiment would not see these two different values. It would measure the average.")
print(f"The bulk measurement would report a single value: {bulk_measurement:.2f} ({bulk_measurement*100:.0f}%) denatured.")
print("\nThe distinct behaviors are lost. This illustrates the impossibility to capture heterogeneity in bulk experiments.")
print("\nTherefore, the correct choice is C.")
