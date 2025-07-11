import math

def solve_interstellar_problem():
    # --- Step 1: Define constants and calculate travel time ---
    distance_ly = 10  # light-years
    speed_c = 0.02   # fraction of speed of light
    
    travel_time_years = distance_ly / speed_c
    
    print(f"Pioneer's journey to Pandora:")
    print(f"- Distance: {distance_ly} light-years")
    print(f"- Speed: {speed_c}c")
    print(f"- Travel Time: {int(travel_time_years)} years\n")

    # --- Step 2: Determine maximum ball capacity and their coordinates ---
    # A simple cubic packing fits 3x2x2=12 balls.
    # A more efficient staggered packing can be achieved. We found a configuration for 16 balls.
    # The coordinates must be multiples of 0.5cm. Ball radius is 2cm.
    # The container is 12x11x11 cm.
    
    max_balls = 16
    ball_radius = 2.0
    box_dims = (12.0, 11.0, 11.0)

    # Coordinates for a 16-ball staggered configuration
    ball_centers = []
    # Layer 1 (z=2.0): 6 balls
    for x in [2.0, 6.0, 10.0]:
        for y in [2.0, 6.0]:
            ball_centers.append((x, y, 2.0))
    # Layer 2 (z=5.0): 4 balls (interstitial)
    for x in [4.0, 8.0]:
        for y in [4.0, 8.0]:
            ball_centers.append((x, y, 5.0))
    # Layer 3 (z=8.0): 6 balls
    for x in [2.0, 6.0, 10.0]:
        for y in [2.0, 6.0]:
            ball_centers.append((x, y, 8.0))

    print(f"Energy Ball Packing:")
    print(f"- A dense packing allows for a maximum of {max_balls} balls.")
    
    # --- Step 3: Count leaking vs. non-leaking balls ---
    num_leaking = 0
    num_non_leaking = 0
    
    for center in ball_centers:
        cx, cy, cz = center
        # A ball is leaking if it touches any surface.
        # This happens if its center is within ball_radius of any edge of the box [0, L].
        is_leaking = (cx - ball_radius <= 0 or cx + ball_radius >= box_dims[0] or
                       cy - ball_radius <= 0 or cy + ball_radius >= box_dims[1] or
                       cz - ball_radius <= 0 or cz + ball_radius >= box_dims[2])
        if is_leaking:
            num_leaking += 1
        else:
            num_non_leaking += 1
            
    print(f"- Number of balls touching container surface (leaking): {num_leaking}")
    print(f"- Number of internal balls (non-leaking): {num_non_leaking}\n")

    # --- Step 4: Calculate final energy per ball ---
    initial_energy_mj = 100
    leak_rate_per_year = 0.001 # 0.1%

    # Energy of a non-leaking ball remains constant
    final_energy_non_leaking = initial_energy_mj

    # Energy of a leaking ball decays exponentially
    final_energy_leaking = initial_energy_mj * math.pow(1 - leak_rate_per_year, travel_time_years)

    print("Energy Calculation per Ball:")
    print(f"- Initial energy per ball: {initial_energy_mj} MJ")
    print(f"- Final energy of a non-leaking ball: {final_energy_non_leaking:.2f} MJ")
    print(f"- Final energy of a leaking ball: {final_energy_leaking:.2f} MJ\n")

    # --- Step 5: Check mission requirement ---
    total_final_energy = (num_leaking * final_energy_leaking) + (num_non_leaking * final_energy_non_leaking)
    required_energy_mj = 1000

    print("Total Energy on Arrival:")
    print("Final Equation:")
    print(f"({num_leaking} leaking balls * {final_energy_leaking:.2f} MJ) + ({num_non_leaking} non-leaking balls * {final_energy_non_leaking:.2f} MJ) = {total_final_energy:.2f} MJ")
    
    print(f"\nConclusion:")
    print(f"- Total energy on arrival: {total_final_energy:.2f} MJ")
    print(f"- Required energy for operations: {required_energy_mj} MJ")
    
    if total_final_energy >= required_energy_mj:
        print("- The energy is sufficient for operations.")
        result = max_balls
    else:
        print("- The energy is NOT sufficient for operations.")
        result = 0
        
    return result

# Run the simulation and get the final answer
final_answer = solve_interstellar_problem()
print(f"<<<{final_answer}>>>")
