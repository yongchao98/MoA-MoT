import math

def solve_cutting_problem():
    """
    Calculates the minimum number of cuts to divide a 4x4x4 cube into 1x1x1 cubes
    with a knife that has a cutting depth of 2cm.
    """
    knife_depth = 2
    total_cuts = 0
    cut_costs = []

    print("--- Optimal Cutting Strategy ---")
    print("The strategy is to reduce the height of the pieces as quickly as possible to minimize the cost of subsequent cuts.")
    print("Let's assume cuts are vertical, so height is the Y-dimension of the pieces.\n")

    # Step 1: Cut the central Y-plane (y=2)
    height_step1 = 4
    cost_step1 = math.ceil(height_step1 / knife_depth)
    total_cuts += cost_step1
    cut_costs.append(cost_step1)
    print(f"1. Cut plane y=2: Piece height is {height_step1}cm. Cost = {cost_step1} cuts.")

    # Subsequent pieces now have a height of 2cm
    height_step2 = 2

    # Step 2: Cut the central X and Z planes (x=2, z=2)
    cost_step2_x = math.ceil(height_step2 / knife_depth)
    total_cuts += cost_step2_x
    cut_costs.append(cost_step2_x)
    print(f"2. Cut plane x=2: Piece height is {height_step2}cm. Cost = {cost_step2_x} cut.")

    cost_step2_z = math.ceil(height_step2 / knife_depth)
    total_cuts += cost_step2_z
    cut_costs.append(cost_step2_z)
    print(f"3. Cut plane z=2: Piece height is {height_step2}cm. Cost = {cost_step2_z} cut.")

    # Step 3: Cut the remaining Y-planes (y=1, y=3)
    # These are 2 planes, each cutting through pieces of height 2cm
    cost_step3 = 2 * math.ceil(height_step2 / knife_depth)
    total_cuts += cost_step3
    cut_costs.append(cost_step3)
    print(f"4. Cut planes y=1 & y=3: Piece height is {height_step2}cm. Cost = {cost_step3} cuts.")

    # Subsequent pieces now have a height of 1cm
    height_step4 = 1

    # Step 4: Cut the remaining X and Z planes
    # There are 2 remaining X-planes (x=1, x=3) and 2 remaining Z-planes (z=1, z=3).
    # All these cuts are through pieces with a height of 1cm.
    cost_step4_x = 2 * math.ceil(height_step4 / knife_depth)
    total_cuts += cost_step4_x
    cut_costs.append(cost_step4_x)
    print(f"5. Cut planes x=1 & x=3: Piece height is {height_step4}cm. Cost = {cost_step4_x} cuts.")
    
    cost_step4_z = 2 * math.ceil(height_step4 / knife_depth)
    total_cuts += cost_step4_z
    cut_costs.append(cost_step4_z)
    print(f"6. Cut planes z=1 & z=3: Piece height is {height_step4}cm. Cost = {cost_step4_z} cuts.")
    
    print("\n--- Final Calculation ---")
    equation_parts = [str(cost) for cost in cut_costs]
    equation = " + ".join(equation_parts)
    print(f"The minimum number of cuts is the sum of costs from each step:")
    print(f"{equation} = {total_cuts}")

solve_cutting_problem()
<<<10>>>