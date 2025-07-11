import math

def solve_horse_area():
    """
    Calculates the area a horse can reach on a rope with taxi-cab length,
    avoiding a house-shaped obstacle.
    """
    # Rope length
    L = 7.0 / 2.0

    # 1. Calculate the area in the three unobstructed quadrants (Q1, Q2, Q4).
    # The total area of a taxi-cab circle of radius L is 2 * L^2.
    # The three quadrants represent 3/4 of this total area.
    area_unobstructed = (3.0 / 4.0) * 2 * L**2
    
    # 2. Calculate the reachable area in the obstructed quadrant (Q3).
    # The rope must pivot around the outer corners of the house.
    # The relevant corners are C1=(-2,0) and C5=(0,-2).

    # Pivot at C1=(-2, 0)
    # Taxi-cab distance from origin (0,0) to C1
    dist_O_C1 = abs(-2 - 0) + abs(0 - 0)
    # Remaining rope length
    r1 = L - dist_O_C1
    # The area reachable from C1 is a taxi-cab circle of radius r1.
    # The part that extends into Q3 (y<0) is half of this circle's area.
    # Area of a full taxi-cab circle is 2 * r^2. Half is r^2.
    area_from_C1 = r1**2

    # Pivot at C5=(0, -2)
    # Taxi-cab distance from origin (0,0) to C5
    dist_O_C5 = abs(0 - 0) + abs(-2 - 0)
    # Remaining rope length
    r5 = L - dist_O_C5
    # The area reachable from C5 is a taxi-cab circle of radius r5.
    # The part that extends into Q3 (x<0) is half of this circle's area.
    area_from_C5 = r5**2
    
    # The regions from C1 and C5 are disjoint, so we sum their areas.
    area_q3 = area_from_C1 + area_from_C5

    # 3. Calculate the total area.
    total_area = area_unobstructed + area_q3

    # 4. Print the explanation and the final equation.
    print("The total reachable area is the sum of the area in the unobstructed quadrants and the area reachable in the obstructed quadrant.")
    print("\nStep 1: Area in unobstructed quadrants (Q1, Q2, Q4)")
    print(f"This is 3/4 of the total area of a taxi-cab circle with radius {L}.")
    print(f"Area_unobstructed = (3/4) * 2 * ({L})^2 = {area_unobstructed}")

    print("\nStep 2: Area in the obstructed quadrant (Q3)")
    print("This area is formed by the rope pivoting around the house's corners.")
    print(f"Pivoting at (-2,0): Remaining rope = {L} - 2 = {r1}. Reachable area in Q3 = ({r1})^2 = {area_from_C1}")
    print(f"Pivoting at (0,-2): Remaining rope = {L} - 2 = {r5}. Reachable area in Q3 = ({r5})^2 = {area_from_C5}")
    print(f"Total area in Q3 = {area_from_C1} + {area_from_C5} = {area_q3}")

    print("\nStep 3: Final Calculation")
    print("The final equation for the total area is:")
    print(f"Total Area = (Area in Q1,Q2,Q4) + (Area from C1) + (Area from C5)")
    print(f"Total Area = (3/4 * 2 * {L}^2) + ({r1}^2) + ({r5}^2)")
    print(f"Total Area = {area_unobstructed} + {area_from_C1} + {area_from_C5}")
    print(f"Total Area = {total_area}")

solve_horse_area()