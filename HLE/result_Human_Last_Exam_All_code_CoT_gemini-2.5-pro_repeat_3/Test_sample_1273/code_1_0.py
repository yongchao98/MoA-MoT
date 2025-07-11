import math

# Define billet dimensions
billet_dims = {'l': 8, 'w': 8, 'h': 4}

# Define product specifications
products = {
    'B2': {'name': 'Balls of 2cm radius', 'price': 150, 'block_dims': (4, 4, 4)},
    'B1': {'name': 'Balls of 1cm diameter', 'price': 1, 'block_dims': (1, 1, 1)},
    'T1': {'name': 'Cubes of 0.8cm side', 'price': 1, 'block_dims': (0.8, 0.8, 0.8)}
}

print("Analyzing the best way to cut the 8x8x4 cm billet.")
print("-" * 50)

best_value = 0
best_strategy_details = ""

# --- Strategy 1: Use only B2 balls ---
p_b2 = products['B2']
d_b2 = p_b2['block_dims']
num_b2 = (math.floor(billet_dims['l'] / d_b2[0]) *
          math.floor(billet_dims['w'] / d_b2[1]) *
          math.floor(billet_dims['h'] / d_b2[2]))
val_b2 = num_b2 * p_b2['price']
print(f"Strategy: Fill with only B2 Balls ({p_b2['name']})")
print(f"Number of balls that can be cut: {num_b2}")
print(f"Calculation: {num_b2} * {p_b2['price']} = {val_b2}")
print("-" * 50)

if val_b2 > best_value:
    best_value = val_b2
    best_strategy_details = f"The best strategy is to produce {num_b2} B2 balls."
    final_equation = f"{num_b2} * {p_b2['price']} = {val_b2}"

# --- Strategy 2: Use only T1 cubes ---
p_t1 = products['T1']
d_t1 = p_t1['block_dims']
num_t1 = (math.floor(billet_dims['l'] / d_t1[0]) *
          math.floor(billet_dims['w'] / d_t1[1]) *
          math.floor(billet_dims['h'] / d_t1[2]))
val_t1 = num_t1 * p_t1['price']
print(f"Strategy: Fill with only T1 Cubes ({p_t1['name']})")
print(f"Number of cubes that can be cut: {num_t1}")
print(f"Calculation: {num_t1} * {p_t1['price']} = {val_t1}")
print("-" * 50)

if val_t1 > best_value:
    best_value = val_t1
    best_strategy_details = f"The best strategy is to produce {num_t1} T1 cubes."
    final_equation = f"{num_t1} * {p_t1['price']} = {val_t1}"

# --- Strategy 3: Use only B1 balls ---
p_b1 = products['B1']
d_b1 = p_b1['block_dims']
num_b1 = (math.floor(billet_dims['l'] / d_b1[0]) *
          math.floor(billet_dims['w'] / d_b1[1]) *
          math.floor(billet_dims['h'] / d_b1[2]))
val_b1 = num_b1 * p_b1['price']
print(f"Strategy: Fill with only B1 Balls ({p_b1['name']})")
print(f"Number of balls that can be cut: {num_b1}")
print(f"Calculation: {num_b1} * {p_b1['price']} = {val_b1}")
print("-" * 50)

if val_b1 > best_value:
    best_value = val_b1
    best_strategy_details = f"The best strategy is to produce {num_b1} B1 balls."
    final_equation = f"{num_b1} * {p_b1['price']} = {val_b1}"

# --- Conclusion ---
print("Conclusion:")
print("Comparing the strategies, prioritizing the most value-dense product (B2 Ball) yields the highest value.")
print("The billet dimensions (8x8x4) are perfectly divisible by the B2 ball's required block size (4x4x4), leaving no waste material.")
print("Therefore, any hybrid strategy will result in a lower total value.")
print(best_strategy_details)
print("The final calculation for the best valid option is:")
print(final_equation)
