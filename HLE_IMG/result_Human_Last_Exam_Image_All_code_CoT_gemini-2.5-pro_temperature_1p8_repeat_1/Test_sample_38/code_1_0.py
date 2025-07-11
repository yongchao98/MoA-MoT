import math

print("This program determines which honeycomb tiling has the lowest Poisson's ratio.")
print("The principle is that a more re-entrant (auxetic) structure has a lower Poisson's ratio.")
print("We model the 're-entrantness' by the ratio of the side lengths b/a.")
print("A higher b/a ratio corresponds to a more re-entrant shape and a lower Poisson's ratio.")
print("\nCalculating the b/a ratio for each option:\n")

options = {
    'A': {'label': '(0, 1)', 'a': 0, 'b': 1},
    'B': {'label': '(1, 4)', 'a': 1, 'b': 4},
    'C': {'label': '(1, sqrt(3))', 'a': 1, 'b': math.sqrt(3)},
    'D': {'label': '(1, 1)', 'a': 1, 'b': 1},
    'E': {'label': '(sqrt(3), 1)', 'a': math.sqrt(3), 'b': 1},
    'F': {'label': '(4, 1)', 'a': 4, 'b': 1},
    'G': {'label': '(1, 0)', 'a': 1, 'b': 0}
}

calculated_ratios = {}
max_ratio = -1.0
best_option = None

# Iterate through options and perform calculations
for key, values in options.items():
    a = values['a']
    b = values['b']
    
    if a == 0:
        ratio = float('inf')
        print(f"Option {key} {values['label']}: ratio = b/a = {b}/{a}. This is infinitely large.")
    else:
        ratio = b / a
        print(f"Option {key} {values['label']}: ratio = b/a = {b:.3f}/{a:.3f} = {ratio:.3f}")
        
    calculated_ratios[key] = ratio
    
    if ratio > max_ratio:
        max_ratio = ratio
        best_option = key

print("\n--- Conclusion ---")
print("Comparing the ratios, the highest value corresponds to the most re-entrant structure.")
print(f"The maximum ratio belongs to Option {best_option}, which is {options[best_option]['label']}.")
print(f"Therefore, the tiling for {options[best_option]['label']} will have the lowest Poisson's ratio.")
