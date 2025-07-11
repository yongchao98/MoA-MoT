import math

def solve_container_design():
    """
    Solves the Pioneer probe container design problem by analyzing box and cylinder
    options to find the one with the minimum total cost.
    """
    # --- Step 1: Define Constants and Initial Calculations ---
    BALL_DIAMETER = 4.0  # cm
    BALL_ENERGY = 30  # MJ
    BALL_COST = 1000  # USD

    CONTAINER_MATERIAL_COST = 200  # USD per cm^2
    MAX_SURFACE_AREA = 1000  # cm^2
    MIN_TOTAL_ENERGY = 1000  # MJ

    # Calculate the minimum number of energy balls required
    min_balls_needed = math.ceil(MIN_TOTAL_ENERGY / BALL_ENERGY)
    cost_balls = min_balls_needed * BALL_COST

    print("--- Design & Cost Analysis ---")
    print(f"1. Minimum energy balls required: ceil({MIN_TOTAL_ENERGY} MJ / {BALL_ENERGY} MJ/ball) = {min_balls_needed} balls")

    # --- Step 2: Analyze Container Options ---
    # The goal is to find a container that holds 34 balls with SA <= 1000 cm^2.

    # Box Analysis: The most compact grid to hold >34 balls is 3x3x4=36 balls.
    # This requires minimum dimensions of L=12, W=12, H=16 cm.
    # SA_box = 2 * (12*12 + 12*16 + 16*12) = 1056 cm^2.
    # Since 1056 > 1000, this is not a valid design. Other grid layouts are less efficient.
    print("\n2. Container Design Choice:")
    print("   - Box Container: Not feasible. The most compact design (for 3x3x4 balls) has a surface area of 1056 cm^2, which exceeds the 1000 cm^2 limit.")
    print("   - Cylinder Container: Feasible. The analysis below shows a valid design.")

    # Cylinder Analysis: Packing balls in layers. We found the optimal packing is 7 balls/layer.
    # To hold 34 balls, we need ceil(34/7) = 5 layers.
    # This design requires a radius R=6.0 cm and a height H=20.0 cm.
    R_optimal = 6.0
    H_optimal = 20.0
    print("\n   The optimal cylinder design found has the following parameters:")
    print(f"   - Packing: 7 balls per layer for 5 layers")
    print(f"   - Container Radius (R): {R_optimal:.1f} cm")
    print(f"   - Container Height (H): {H_optimal:.1f} cm")

    # --- Step 3: Calculate Final Cost for the Optimal Design ---
    # Surface area of the optimal cylinder
    sa_cylinder = 2 * math.pi * R_optimal**2 + 2 * math.pi * R_optimal * H_optimal

    # Cost of the container
    cost_container = sa_cylinder * CONTAINER_MATERIAL_COST
    
    # Total cost
    total_cost = cost_container + cost_balls

    # --- Step 4: Output the results and final equation ---
    print("\n3. Final Cost Calculation:")
    print(f"   - Surface Area (SA) = 2 * pi * R^2 + 2 * pi * R * H")
    final_sa_equation = f"   - SA = 2 * pi * {R_optimal:.1f}^2 + 2 * pi * {R_optimal:.1f} * {H_optimal:.1f} = {sa_cylinder:.4f} cm^2"
    print(final_sa_equation)
    
    print(f"\n   - Container Cost = Surface Area * Cost per cm^2")
    final_container_cost_equation = f"   - Container Cost = {sa_cylinder:.4f} * ${CONTAINER_MATERIAL_COST} = ${cost_container:.2f}"
    print(final_container_cost_equation)

    print(f"\n   - Energy Balls Cost = Number of Balls * Cost per Ball")
    final_ball_cost_equation = f"   - Energy Balls Cost = {min_balls_needed} * ${BALL_COST} = ${cost_balls:.2f}"
    print(final_ball_cost_equation)

    print("\n4. Final Total Cost (C):")
    print(f"   C = (Container Cost) + (Energy Balls Cost)")
    print(f"   C = ${cost_container:.2f} + ${cost_balls:.2f}")
    final_total_cost_equation = f"   C = ${total_cost:.2f}"
    print(final_total_cost_equation)
    
    # Final answer as requested by the format
    final_answer = round(total_cost)
    print(f"\n<<<C={final_answer}>>>")

solve_container_design()