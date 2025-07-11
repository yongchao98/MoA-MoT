# --- Experimental Parameters ---

# Number of wells where paired heavy and light chains were successfully amplified
successful_wells = 4

# Total number of wells that received a sorted doublet
total_sorted_wells = 95


# --- Calculation ---

# The efficiency is the ratio of successful outcomes to the total number of attempts.
efficiency = successful_wells / total_sorted_wells


# --- Output ---

print("Experimental Efficiency Calculation:")
print(f"Successful wells: {successful_wells}")
print(f"Total sorted wells: {total_sorted_wells}")
print("\nEquation:")
print(f"{successful_wells} / {total_sorted_wells} = {efficiency}")
print(f"\nThe experimental efficiency is {efficiency:.4f}, or {efficiency:.2%}.")
print("\nThis low efficiency suggests that the event of interest is very rare, or the detection method captures many non-specific events.")