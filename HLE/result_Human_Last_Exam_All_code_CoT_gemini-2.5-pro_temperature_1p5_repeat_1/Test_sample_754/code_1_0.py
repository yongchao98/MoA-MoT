import math

def solve_interstellar_problem():
    # --- Step 1: Define constants and calculate travel time ---
    distance_ly = 10.0  # light-years
    speed_c = 0.02     # as a fraction of the speed of light c
    required_energy_mj = 1000.0
    
    initial_energy_per_ball_mj = 100.0
    leak_rate_per_year = 0.001  # 0.1%

    box_dims = (12, 11, 11) # cm
    ball_radius = 2.0 # cm
    
    travel_time_years = distance_ly / speed_c
    
    print("--- Problem Analysis ---")
    print(f"The probe needs to travel {distance_ly} light-years at {speed_c}c.")
    print(f"Equation: Travel Time = Distance / Speed = {distance_ly} / {speed_c}")
    print(f"Calculated travel time: {travel_time_years:.0f} years.\n")

    # --- Step 2: Calculate energy remaining per ball ---
    energy_factor = (1 - leak_rate_per_year) ** travel_time_years
    energy_leaking_ball = initial_energy_per_ball_mj * energy_factor
    energy_non_leaking_ball = initial_energy_per_ball_mj
    
    print("--- Energy Calculation per Ball ---")
    print(f"Energy of a non-leaking ball remains at {energy_non_leaking_ball:.2f} MJ.")
    print(f"Energy of a leaking ball decays over {travel_time_years:.0f} years.")
    print(f"Equation: Final Energy = Initial Energy * (1 - Leak Rate) ^ Time = {initial_energy_per_ball_mj:.0f} * (1 - {leak_rate_per_year}) ^ {travel_time_years:.0f}")
    print(f"Final energy of a leaking ball: {energy_leaking_ball:.2f} MJ.\n")
    
    # --- Step 3: Analyze Packing Strategies ---
    print("--- Sphere Packing Analysis ---")
    print("The container is 12x11x11 cm. Ball radius is 2 cm.")
    print("A ball's center (x,y,z) must be in the range [2, 10] for x, [2, 9] for y, and [2, 9] for z.")
    
    # Strategy A: Simple Cubic Packing (Greedy placement starting from the corner)
    # Ball centers are placed on a 4cm grid.
    num_x_a = math.floor((box_dims[0] - 2 * ball_radius) / (2 * ball_radius)) + 1 # (12-4)/4 + 1 = 3
    num_y_a = math.floor((box_dims[1] - 2 * ball_radius) / (2 * ball_radius)) + 1 # (11-4)/4 + 1 = 2
    num_z_a = math.floor((box_dims[2] - 2 * ball_radius) / (2 * ball_radius)) + 1 # (11-4)/4 + 1 = 2
    total_balls_a = num_x_a * num_y_a * num_z_a
    # Internal balls are those not touching any wall. In this packing, only the one at center (6,6,6) is internal.
    internal_balls_a = 1
    leaking_balls_a = total_balls_a - internal_balls_a
    total_energy_a = (internal_balls_a * energy_non_leaking_ball) + (leaking_balls_a * energy_leaking_ball)

    print("\nStrategy A: Simple Cubic Packing")
    print(f"This packing fits {total_balls_a} balls: {internal_balls_a} non-leaking and {leaking_balls_a} leaking.")
    print(f"Total Energy = {internal_balls_a} * {energy_non_leaking_ball:.2f} MJ + {leaking_balls_a} * {energy_leaking_ball:.2f} MJ")
    print(f"Total energy from Strategy A: {total_energy_a:.2f} MJ.")
    
    # Strategy B: Maximize Internal Balls
    # Internal center space: (2, 10) x (2, 9) x (2, 9). This is a volume of (10-2)x(9-2)x(9-2) = 8x7x7. No, its (9.5-2.5)x(8.5-2.5)x(8.5-2.5)=7x6x6
    # Size of internal space for centers: x: 7cm, y: 6cm, z: 6cm
    num_x_b = math.floor(7.0 / (2 * ball_radius)) + 1 # floor(7/4)+1=2
    num_y_b = math.floor(6.0 / (2 * ball_radius)) + 1 # floor(6/4)+1=2
    num_z_b = math.floor(6.0 / (2 * ball_radius)) + 1 # floor(6/4)+1=2
    total_balls_b = num_x_b * num_y_b * num_z_b
    internal_balls_b = total_balls_b
    leaking_balls_b = 0
    total_energy_b = internal_balls_b * energy_non_leaking_ball

    print("\nStrategy B: Maximize Non-Leaking (Internal) Balls")
    print("This strategy focuses on packing balls in the central volume away from walls.")
    print(f"This packing fits {total_balls_b} balls: {internal_balls_b} non-leaking and {leaking_balls_b} leaking.")
    print(f"Total Energy = {internal_balls_b} * {energy_non_leaking_ball:.2f} MJ + {leaking_balls_b} * {energy_leaking_ball:.2f} MJ")
    print(f"Total energy from Strategy B: {total_energy_b:.2f} MJ.")
    
    # --- Step 4 & 5: Conclusion and Final Answer ---
    print("\n--- Conclusion ---")
    max_energy = max(total_energy_a, total_energy_b)
    
    if max_energy == total_energy_a:
        max_balls = total_balls_a
    else:
        max_balls = total_balls_b
        
    print(f"The maximum possible stored energy upon arrival is {max_energy:.2f} MJ.")
    print(f"The required energy for operations is {required_energy_mj:.2f} MJ.")
    
    if max_energy >= required_energy_mj:
        print("The container IS sufficient.")
        final_answer = max_balls
    else:
        print("The container IS NOT sufficient.")
        final_answer = 0
        
    print(f"\nFinal Answer: {final_answer}")
    return final_answer

solve_interstellar_problem()