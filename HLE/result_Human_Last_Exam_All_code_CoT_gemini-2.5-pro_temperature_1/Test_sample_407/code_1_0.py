import math

# This script calculates the selection coefficient 's', a measure of fitness cost,
# from a hypothetical competition experiment between a hybrid and a parental yeast strain.

# --- Hypothetical Experimental Data ---

# R0 is the initial ratio of hybrid cells to parental cells.
# Let's assume we start with an equal number of both.
initial_ratio_R0 = 1.0 

# Rt is the final ratio of hybrid to parental cells after competition.
# If the hybrid is less fit, its proportion will decrease.
# Let's say the final ratio is 0.5 (one hybrid cell for every two parental cells).
final_ratio_Rt = 0.5

# t is the number of generations the competition ran for.
generations_t = 10

# --- Calculation ---
# The selection coefficient 's' can be calculated using the formula:
# s = [ln(Rt / R0)] / t
# where 'ln' is the natural logarithm.

# A negative 's' indicates a fitness cost for the hybrid.
selection_coefficient_s = (math.log(final_ratio_Rt / initial_ratio_R0)) / generations_t

# --- Output the results ---
print("Calculating the selection coefficient 's' to measure the cost of gene flow.")
print("-" * 60)
print(f"Initial Ratio (Hybrid/Parental), R0 = {initial_ratio_R0}")
print(f"Final Ratio (Hybrid/Parental), Rt   = {final_ratio_Rt}")
print(f"Number of Generations, t           = {generations_t}")
print("-" * 60)
print("Formula: s = [ln(Rt / R0)] / t")
print(f"Calculation: s = [ln({final_ratio_Rt} / {initial_ratio_R0})] / {generations_t}")
print(f"The selection coefficient 's' is: {selection_coefficient_s:.4f}")

if selection_coefficient_s < 0:
    print("\nSince 's' is negative, this indicates a fitness cost for the hybrid strain,")
    print("which is a measured cost due to gene flow.")
else:
    print("\nSince 's' is not negative, no fitness cost was observed in this experiment.")
