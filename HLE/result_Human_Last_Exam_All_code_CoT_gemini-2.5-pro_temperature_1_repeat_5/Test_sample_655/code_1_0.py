import math

# --- 1. Define the container and its properties ---
container_type = "box"
L, W, H = 12.0, 12.0, 15.5

# --- 2. Verify the surface area constraint ---
surface_area = 2 * (L*W + L*H + W*H)
is_valid_area = surface_area <= 1050

print(f"Chosen container: {container_type} {L}x{W}x{H}")
print(f"Calculated surface area: {surface_area:.2f} cm^2")
print(f"Surface area is within the 1050 cm^2 limit: {is_valid_area}")
print("-" * 20)

# --- 3. Calculate the number of packed balls based on the plan ---

# Number of large balls (n2)
# A 3x3x3 arrangement of large balls in 4x4x4 cells.
n2_x = math.floor(L / 4)
n2_y = math.floor(W / 4)
n2_z = math.floor(12.0 / 4) # They fit in a 12x12x12 subspace
n2 = n2_x * n2_y * n2_z

print("Packing large balls (radius 2cm, diameter 4cm):")
print(f"Arrangement: {n2_x} x {n2_y} x {n2_z}")
print(f"Number of 2-cm balls (b): {n2}")
print("-" * 20)

# Number of small balls (n1)
# They are packed in the remaining space: a 12x12x3.5 slab.
rem_L, rem_W, rem_H = L, W, H - 12.0
n1_x = math.floor(rem_L / 2)
n1_y = math.floor(rem_W / 2)
n1_z = math.floor(rem_H / 2)
n1 = n1_x * n1_y * n1_z

print("Packing small balls (radius 1cm, diameter 2cm) in remaining space:")
print(f"Remaining space dimensions: {rem_L}x{rem_W}x{rem_H}")
print(f"Arrangement: {n1_x} x {n1_y} x {n1_z}")
print(f"Number of 1-cm balls (a): {n1}")
print("-" * 20)

# --- 4. Calculate the total potential energy ---
energy_n1 = 1  # MJ per 1-cm ball
energy_n2 = 10 # MJ per 2-cm ball
total_energy = (n1 * energy_n1) + (n2 * energy_n2)

print("Calculating total energy:")
print(f"Final equation: ({n1} * {energy_n1} MJ) + ({n2} * {energy_n2} MJ) = {total_energy} MJ")
print("-" * 20)

# --- 5. Format the final answer ---
final_answer_str = f"[{container_type} {L}x{W}x{H}]{n1};{n2}"
print(f"Final Answer String: {final_answer_str}")

# Output the final answer in the required format
print("<<<[box 12.0x12.0x15.5]36;27>>>")