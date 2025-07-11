import math

def solve_puzzle():
    """
    Solves the AGI QA puzzle by modeling the geometry and using the provided data.
    """

    # --- Step 1: Analyze the Data ---
    monte_carlo_points = [730, 740, 735, 732, 739]
    total_points = 1000
    average_points_in_circles = sum(monte_carlo_points) / len(monte_carlo_points)
    
    # Ratio of the area of 9 circles to the total area
    ratio_circles_area = average_points_in_circles / total_points
    
    # Ratio of one green rectangle's area to the total area
    ratio_green_rect_area = 0.04  # "around 4%"

    print("Step 1: Processing the given information")
    print(f"Average points in circles: {average_points_in_circles}")
    print(f"Estimated ratio of circles' area to total area: {ratio_circles_area:.4f}")
    print(f"Given ratio of a green rectangle's area to total area: ~{ratio_green_rect_area:.2f}")
    print("-" * 30)

    # --- Step 2 & 3: Develop and Verify a Geometric Model ---
    # From visual inspection of the image, we establish a geometric model.
    # - Let R be the radius of the yellow circles.
    # - The height of a green rectangle (gh) appears to be equal to the diameter of a circle: gh = 2R
    # - The circles are in a staggered packing. The most stable arrangement gives:
    #   - Total height y = 2R * (1 + sqrt(3))
    #   - Total width x = 7R
    # We will now test the consistency of this model.

    print("Step 2: Proposing and testing a geometric model")
    print("Model based on visual layout:")
    print("  - Outer rectangle width, x = 7 * R")
    print("  - Outer rectangle height, y = 2 * R * (1 + sqrt(3))")
    print("  - Green rectangle height, gh = 2 * R")

    # Let's verify the circles' area ratio with this model.
    # Theoretical ratio = (Area of 9 circles) / (Area of rectangle)
    #                   = (9 * pi * R^2) / (x * y)
    #                   = (9 * pi * R^2) / (7R * 2R * (1 + sqrt(3)))
    #                   = (9 * pi) / (14 * (1 + sqrt(3)))
    
    model_ratio_circles_area = (9 * math.pi) / (14 * (1 + math.sqrt(3)))
    print(f"\nModel's predicted ratio for circles' area: {model_ratio_circles_area:.4f}")
    print(f"This is very close to the Monte Carlo result of {ratio_circles_area:.4f}, so the model is strong.")

    # Now, let's use the model to find the relationship between gw and R.
    # Ratio = (Area of green rectangle) / (Area of rectangle)
    # 0.04  = (gw * gh) / (x * y)
    # 0.04  = (gw * 2R) / (7R * 2R * (1 + sqrt(3)))
    # 0.04  = gw / (7R * (1 + sqrt(3)))
    # Solving for gw/R:
    # gw/R = 0.04 * 7 * (1 + sqrt(3))

    gw_over_R_ratio = ratio_green_rect_area * 7 * (1 + math.sqrt(3))
    print(f"\nFrom the 4% area rule, the model predicts gw/R should be ~{gw_over_R_ratio:.4f}")
    print("-" * 30)

    # --- Step 4: Solve for Integer Dimensions ---
    # We are told R and gw are integers. We need to find an integer pair (R, gw)
    # such that gw/R is approximately 0.765. The simplest rational approximation is 3/4 = 0.75.
    
    print("Step 3: Finding the integer dimensions")
    print("The ratio gw/R is ~0.765. The closest simple fraction is 3/4 = 0.75.")
    print("For gw and R to be integers, R must be a multiple of 4.")
    print("The simplest possible integer solution is R = 4, which gives gw = 3.")

    R = 4
    gw = 3
    gh = 2 * R
    
    print(f"\nChosen values: R = {R} cm, gw = {gw} cm, gh = {gh} cm.")
    print("-" * 30)

    # --- Step 5: Calculate Final Answer ---
    print("Step 4: Calculating the final dimensions of the outer rectangle")
    x = 7 * R
    y = 2 * R * (1 + math.sqrt(3))

    print(f"Calculated length x = 7 * {R} = {x}")
    print(f"Calculated width y = 2 * {R} * (1 + sqrt(3)) = {y:.3f}")

    # Round to the nearest cm as requested
    x_rounded = round(x)
    y_rounded = round(y)
    
    print("\nRounding to the nearest cm:")
    print(f"  Final length x = {x_rounded} cm")
    print(f"  Final width y = {y_rounded} cm")
    
    print("\nThe final answer is...")
    print(f"{x_rounded}:{y_rounded}")

solve_puzzle()