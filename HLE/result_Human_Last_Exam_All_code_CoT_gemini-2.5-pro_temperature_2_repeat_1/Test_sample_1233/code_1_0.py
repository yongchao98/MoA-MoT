import math

def simulate_bacterial_growth():
    """
    Simulates the growth of a slow pathogen vs. a fast contaminant
    to demonstrate why longer incubation is necessary.
    """
    # --- Parameters ---
    # Assume the contaminant grows much faster than the target pathogen.
    # Doubling time in hours (lower is faster)
    contaminant_doubling_time = 2.0
    pathogen_doubling_time = 4.0

    # Initial number of bacteria plated
    initial_count = 100

    # Colony becomes visible when count exceeds this threshold
    visibility_threshold = 1_000_000

    print("Step 1: The lab incubates the plate for 48 hours (2 days).\n")
    time_1 = 48 # hours

    # --- Calculations at 48 hours ---
    # Growth formula: final_count = initial_count * 2^(time / doubling_time)
    contaminant_count_48 = initial_count * (2**(time_1 / contaminant_doubling_time))
    pathogen_count_48 = initial_count * (2**(time_1 / pathogen_doubling_time))
    
    # --- Report results at 48 hours ---
    print(f"--- Analysis at {time_1} hours ---")
    print(f"Contaminant count = {initial_count} * 2^({time_1}/{contaminant_doubling_time}) = {contaminant_count_48:,.0f}")
    print(f"Pathogen count = {initial_count} * 2^({time_1}/{pathogen_doubling_time}) = {pathogen_count_48:,.0f}")
    
    contaminant_visible_48 = "YES" if contaminant_count_48 > visibility_threshold else "NO"
    pathogen_visible_48 = "YES" if pathogen_count_48 > visibility_threshold else "NO"

    print(f"\nIs Contaminant visible (> {visibility_threshold:,.0f})? {contaminant_visible_48}")
    print(f"Is Pathogen visible (> {visibility_threshold:,.0f})? {pathogen_visible_48}")
    print("\nResult at 48 hours: The technician sees the fast-growing contaminant but not the pathogen. This leads to a presumptive negative result.")
    
    print("\n----------------------------------------------------------\n")
    
    print("Step 2: A potential solution is to incubate the existing plates for longer.\n")
    time_2 = 72 # hours

    # --- Calculations at 72 hours ---
    contaminant_count_72 = initial_count * (2**(time_2 / contaminant_doubling_time))
    pathogen_count_72 = initial_count * (2**(time_2 / pathogen_doubling_time))

    # --- Report results at 72 hours ---
    print(f"--- Analysis at {time_2} hours ---")
    print(f"Contaminant count = {initial_count} * 2^({time_2}/{contaminant_doubling_time}) = {contaminant_count_72:,.0f}")
    print(f"Pathogen count = {initial_count} * 2^({time_2}/{pathogen_doubling_time}) = {pathogen_count_72:,.0f}")
    
    contaminant_visible_72 = "YES" if contaminant_count_72 > visibility_threshold else "NO"
    pathogen_visible_72 = "YES" if pathogen_count_72 > visibility_threshold else "NO"
    
    print(f"\nIs Contaminant visible (> {visibility_threshold:,.0f})? {contaminant_visible_72}")
    print(f"Is Pathogen visible (> {visibility_threshold:,.0f})? {pathogen_visible_72}")

    print("\nResult at 72 hours: The slower-growing pathogen has now reached a visible colony size.")
    print("\nConclusion: By incubating the plates for longer, the first lab could have potentially recovered the pathogen despite the overgrowth.")


if __name__ == '__main__':
    simulate_bacterial_growth()
