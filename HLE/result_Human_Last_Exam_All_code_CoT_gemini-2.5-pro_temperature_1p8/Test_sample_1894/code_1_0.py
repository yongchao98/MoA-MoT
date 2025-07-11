import math

def solve_visibility():
    """
    Analyzes the visibility of 5 colored balls in a room from Joe's perspective in the doorway.
    This function will check each ball one by one, showing the calculations.
    """
    
    # --- Room and Object Setup ---
    # Joe's viewing area: x in [4.5, 7.5], y in [0, 0.25], z_eye = 5.0
    # A ball is a sphere with a 0.5ft diameter, so its radius is 0.25ft.

    balls = {
        "Red": {"pos": (11.75, 0.25, 0.25), "visible": False, "reason": ""},
        "Blue": {"pos": (11.75, 11.75, 2.75), "visible": False, "reason": ""},
        "Yellow": {"pos": (0.25, 11.75, 0.25), "visible": False, "reason": ""},
        "Green": {"pos": (0.25, 0.25, 7.25), "visible": False, "reason": ""},
        "Purple": {"pos": (11.75, 4.25, 0.25), "visible": False, "reason": ""}
    }
    
    # --- Analysis for each ball ---

    print("Analyzing visibility for each ball:\n")

    # 1. Red Ball (SE corner, on floor)
    # Obstruction: The door, hinged at (7.5, 0), open 90 degrees inward.
    # The door is a vertical plane at x=7.5, for y between 0 and 3.
    # Joe must stand to the left of the hinge (vx < 7.5) to see past it.
    # Let's take Joe's best position: as far left as possible.
    viewer_pos_red = (4.5, 0.25, 5.0)
    target_pos_red = balls["Red"]["pos"]
    print("--- Checking Red Ball ---")
    print(f"Joe's position: {viewer_pos_red}, Red Ball position: {target_pos_red}")
    print("The primary obstruction is the open door, a plane at x = 7.5.")
    # The line of sight is from vx=4.5 to tx=11.75. It must cross x=7.5.
    # The line's y-coordinate is almost constant at y=0.25.
    # At x = 7.5, the line's y-coordinate will be extremely close to 0.25.
    # The door extends from y=0 to y=3.
    # Since the line of sight at the door plane (x=7.5) has a y-value of ~0.25,
    # which is within the door's y-range [0, 3], the view is blocked.
    balls["Red"]["reason"] = "Blocked by the open door at x=7.5."
    print(f"Result: The Red ball is NOT visible. {balls['Red']['reason']}\n")


    # 2. Green Ball (SW corner, on top of 7ft tall bookshelf)
    # Obstruction: The bookshelf itself. Bookshelf is 0<=x<=1, 0<=y<=4, 0<=z<=7.
    # The ball is above Joe's eye level (7.25ft vs 5ft), so he looks up.
    # Joe's best chance is to be as far as possible to see over the shelf.
    viewer_pos_green = (7.5, 0.25, 5.0) 
    target_pos_green = balls["Green"]["pos"]
    print("--- Checking Green Ball ---")
    print(f"Joe's position: {viewer_pos_green}, Green Ball position: {target_pos_green}")
    print("The primary obstruction is the bookshelf, which has a top edge at z = 7.")
    print("We check if the line of sight passes over the front edge of the shelf (x=1).")
    
    vx, _, vz = viewer_pos_green
    tx, _, tz = target_pos_green
    # Equation for t: vx + t * (tx - vx) = 1
    t = (1.0 - vx) / (tx - vx)
    print(f"Parameter t for intersection with plane x=1: (1.0 - {vx}) / ({tx} - {vx}) = {-6.5} / {-7.25} = {t:.4f}")
    
    # Z coordinate at intersection: z = vz + t * (tz - vz)
    intersect_z = vz + t * (tz - vz)
    print(f"Z-coordinate at x=1: {vz} + {t:.4f} * ({tz} - {vz}) = {vz} + {t:.4f} * {2.25} = {intersect_z:.4f}")

    if intersect_z > 7.0:
        balls["Green"]["visible"] = True
        balls["Green"]["reason"] = f"Line of sight passes at z={intersect_z:.4f}, which is above the shelf's top edge at z=7.0."
    else:
        balls["Green"]["reason"] = f"Line of sight passes at z={intersect_z:.4f}, which is below the shelf's top edge at z=7.0."
    print(f"Result: The Green ball IS visible. {balls['Green']['reason']}\n")

    # 3. Yellow Ball (NW corner, on floor)
    # Obstruction: The bookshelf. 0<=x<=1, 0<=y<=4.
    # Joe's best position is far to the right.
    viewer_pos_yellow = (7.5, 0.25, 5.0)
    target_pos_yellow = balls["Yellow"]["pos"]
    print("--- Checking Yellow Ball ---")
    print(f"Joe's position: {viewer_pos_yellow}, Yellow Ball position: {target_pos_yellow}")
    print("The primary obstruction is the bookshelf body (from y=0 to y=4).")
    print("We check where the line of sight crosses the bookshelf's front plane (x=1).")
    vx, vy, _ = viewer_pos_yellow
    tx, ty, _ = target_pos_yellow
    # t is the same as for the Green Ball since x coordinates are the same
    t = (1.0 - vx) / (tx - vx)
    print(f"Parameter t for intersection with plane x=1: (1.0 - {vx}) / ({tx} - {vx}) = {-6.5} / {-7.25} = {t:.4f}")
    # Y coordinate at intersection: y = vy + t * (ty - vy)
    intersect_y = vy + t * (ty - vy)
    print(f"Y-coordinate at x=1: {vy} + {t:.4f} * ({ty} - {vy}) = {vy} + {t:.4f} * {11.5} = {intersect_y:.4f}")

    if intersect_y > 4.0:
        balls["Yellow"]["visible"] = True
        balls["Yellow"]["reason"] = f"Line of sight passes at y={intersect_y:.4f}, which is beyond the bookshelf's side edge at y=4.0."
    else:
        balls["Yellow"]["reason"] = f"Line of sight passes at y={intersect_y:.4f}, blocked by the bookshelf."
    print(f"Result: The Yellow ball IS visible. {balls['Yellow']['reason']}\n")

    # 4. Blue Ball (NE corner, on a table)
    # Obstruction: The wardrobe. Body at 9.5<=x<=12, 4<=y<=8. Top-left corner is at (9.5, 8).
    # Joe's best chance is from the far left of the door.
    viewer_pos_blue = (4.5, 0.25, 5.0)
    target_pos_blue = balls["Blue"]["pos"]
    print("--- Checking Blue Ball ---")
    print(f"Joe's position: {viewer_pos_blue}, Blue Ball position: {target_pos_blue}")
    print("The primary obstruction is the top corner of the wardrobe (y=8) at its front plane (x=9.5).")
    vx, vy, _ = viewer_pos_blue
    tx, ty, _ = target_pos_blue
    # Equation for t: vx + t * (tx - vx) = 9.5
    t = (9.5 - vx) / (tx - vx)
    print(f"Parameter t for intersection with plane x=9.5: (9.5 - {vx}) / ({tx} - {vx}) = {5.0} / {7.25} = {t:.4f}")
    # Y coordinate at intersection: y = vy + t * (ty - vy)
    intersect_y = vy + t * (ty - vy)
    print(f"Y-coordinate at x=9.5: {vy} + {t:.4f} * ({ty} - {vy}) = {vy} + {t:.4f} * {11.5} = {intersect_y:.4f}")
    
    if intersect_y > 8.0:
        balls["Blue"]["visible"] = True
        balls["Blue"]["reason"] = f"Line of sight passes at y={intersect_y:.4f}, which is just above the wardrobe's top edge at y=8.0."
    else:
        balls["Blue"]["reason"] = f"Line of sight passes at y={intersect_y:.4f}, blocked by the wardrobe."
    print(f"Result: The Blue ball IS visible. {balls['Blue']['reason']}\n")

    # 5. Purple Ball (In SE corner of the wardrobe)
    # Obstruction: The wardrobe's own front panel. Opening is from y=4 to y=8.
    # Joe's best chance is from far left.
    viewer_pos_purple = (4.5, 0.25, 5.0)
    target_pos_purple = balls["Purple"]["pos"]
    print("--- Checking Purple Ball ---")
    print(f"Joe's position: {viewer_pos_purple}, Purple Ball position: {target_pos_purple}")
    print("The primary obstruction is the wardrobe's front panel below the doors. The opening starts at y=4.")
    vx, vy, _ = viewer_pos_purple
    tx, ty, _ = target_pos_purple
    # t is the same as for the Blue Ball as we check the same plane (x=9.5) from the same viewpoint
    t = (9.5 - vx) / (tx - vx)
    print(f"Parameter t for intersection with plane x=9.5: (9.5 - {vx}) / ({tx} - {vx}) = {5.0} / {7.25} = {t:.4f}")
    # Y coordinate at intersection: y = vy + t * (ty - vy)
    intersect_y = vy + t * (ty - vy)
    print(f"Y-coordinate at x=9.5: {vy} + {t:.4f} * ({ty} - {vy}) = {vy} + {t:.4f} * {4.0} = {intersect_y:.4f}")

    if intersect_y >= 4.0:
        balls["Purple"]["visible"] = True
        balls["Purple"]["reason"] = f"Line of sight passes at y={intersect_y:.4f}, which is within the wardrobe opening (y>=4)."
    else:
        balls["Purple"]["reason"] = f"Line of sight passes at y={intersect_y:.4f}, which is below the opening, hitting the solid front panel."
    print(f"Result: The Purple ball is NOT visible. {balls['Purple']['reason']}\n")
    
    # --- Final Conclusion ---
    visible_balls = [name for name, data in balls.items() if data["visible"]]
    print("--- Conclusion ---")
    print("The balls Joe can see are:")
    for ball in visible_balls:
        print(f"- {ball}")
    
    return visible_balls

visible_balls = solve_visibility()
# The final answer format is specific
final_answer = ", ".join(sorted(visible_balls))
# For example: Blue, Green, Yellow
# The brackets will be used to delimit the final answer.
# This print will be captured by the calling system.
print(f"<<<{final_answer}>>>")