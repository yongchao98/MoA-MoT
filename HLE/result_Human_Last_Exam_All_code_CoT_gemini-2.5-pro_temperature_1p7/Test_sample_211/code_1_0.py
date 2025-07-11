import numpy as np

def calculate_and_print_slope(y1, y2, x1, x2, state_label):
    """Calculates and prints the slope of a line between two points."""
    # The slope 'm' of a line is (change in y) / (change in x)
    slope = (y2 - y1) / (x2 - x1)
    
    print(f"--- {state_label} State ---")
    print(f"The equation for the slope is: m = (diversity_2 - diversity_1) / (latitude_2 - latitude_1)")
    print(f"Using the extreme sites, Tena (Site 1) and Ottawa (Site 2):")
    print(f"  - Diversity at Site 1: {y1}")
    print(f"  - Diversity at Site 2: {y2}")
    print(f"  - Latitude at Site 1: {x1}")
    print(f"  - Latitude at Site 2: {x2}")
    # The prompt requests that each number in the final equation be outputted.
    print(f"Calculation: m = ({y2} - {y1}) / ({x2} - ({x1}))")
    print(f"Resulting {state_label.lower()} slope: {slope:.4f}\n")
    return slope

# --- Site Data ---
# Approximate latitudes for the most equatorial and most polar sites.
lat_tena_ecuador = -0.99  # Latitude for Tena, Ecuador
lat_ottawa_canada = 45.42 # Latitude for Ottawa, Canada

# --- Initial State: A steep latitudinal diversity gradient exists ---
# Let's assume high species richness (alpha diversity) in the tropics
# and low richness in the temperate zone.
initial_diversity_tena = 150 
initial_diversity_ottawa = 25

# --- Ultimate State: Biotic homogenization due to invasive species ---
# Under Hubbell's theory, the invasive out-competes natives via stochastic dominance,
# leading to low diversity everywhere.
ultimate_diversity_tena = 3
ultimate_diversity_ottawa = 2

# Calculate and display the initial slope
initial_slope = calculate_and_print_slope(
    initial_diversity_tena,
    initial_diversity_ottawa,
    lat_tena_ecuador,
    lat_ottawa_canada,
    "Initial"
)

# Calculate and display the ultimate slope
ultimate_slope = calculate_and_print_slope(
    ultimate_diversity_tena,
    ultimate_diversity_ottawa,
    lat_tena_ecuador,
    lat_ottawa_canada,
    "Ultimate"
)

print("--- Conclusion ---")
print("The steep initial slope reflects the strong latitudinal diversity gradient.")
print("After the invasion leads to biotic homogenization, diversity becomes uniformly low.")
print("As a result, the final slope is much flatter, having decreased significantly and approached zero.")
