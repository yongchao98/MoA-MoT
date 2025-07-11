import math

def calculate_k(sides):
    """
    Calculates the number of squares crossed by a triangle with given side projections.
    - sides is a list of tuples, where each tuple is (dx, dy) for a side.
    - The calculation assumes vertices are in distinct squares, leading to 3 overlaps.
    """
    
    total_k = 0
    print("Calculating squares crossed for each side:")
    
    # Calculate k for each side
    side_names = ["AB (leg)", "BC (leg)", "AC (hypotenuse)"]
    k_values = []
    for i, (name, (dx, dy)) in enumerate(zip(side_names, sides)):
        # Number of crossed vertical grid lines is floor(dx)
        # Number of crossed horizontal grid lines is floor(dy)
        v_crossings = math.floor(dx)
        h_crossings = math.floor(dy)
        # k = 1 (start square) + vertical crossings + horizontal crossings
        k = 1 + v_crossings + h_crossings
        k_values.append(k)
        print(f"- Side {name} with |Δx|={dx:.2f}, |Δy|={dy:.2f}:")
        print(f"  k({name}) = 1 + floor({dx:.2f}) + floor({dy:.2f}) = 1 + {v_crossings} + {h_crossings} = {k}")
        total_k += k

    # Account for overlaps at vertices
    # Assuming the 3 vertices lie in 3 different squares, each of these squares is counted twice.
    overlaps = 3
    print("\nTotal sum of squares (with double counting) =", " + ".join(map(str, k_values)), "=", total_k)
    print(f"Subtracting {overlaps} for overlapping squares at the vertices.")
    
    final_k = total_k - overlaps
    
    # The final equation is the sum of k for each side minus the overlaps
    print("\nFinal equation for total distinct squares:")
    # Build the final equation string with all numbers
    equation_parts = []
    for k_val in k_values:
      equation_parts.append(str(k_val))
    
    print(f"k = ({' + '.join(equation_parts)}) - {overlaps} = {final_k}")
    
    return final_k

# Optimal projections based on maximizing floor(dx)+floor(dy) under dx^2+dy^2=L^2
# For L=18 (legs), max(floor(dx)+floor(dy)) = 25, achieved with dx=11.3, dy~14.01 or dx~14.01, dy=11.3
# Using the tan(theta)=1/2 analysis for a more concrete setup
dx_leg1 = 36 / math.sqrt(5)
dy_leg1 = 18 / math.sqrt(5)
dx_leg2 = 18 / math.sqrt(5)
dy_leg2 = 36 / math.sqrt(5)

# The hypotenuse is the vector sum, but since we are maximizing, we align it optimally too.
# For L=18*sqrt(2), dx^2+dy^2=648. Let's find projections for one of the leg configurations.
# if AB = (dx_leg1, dy_leg1) and BC = (dx_leg2, -dy_leg1), AC would be (dx_leg1+dx_leg2, dy_leg1-dy_leg2) NO
# Vector AC = AB+BC, AB perp BC
# AB = (36/sqrt(5), 18/sqrt(5)), BC = (-18/sqrt(5), 36/sqrt(5))
# AC = AB+BC = (18/sqrt(5), 54/sqrt(5))
dx_hyp = 54 / math.sqrt(5)
dy_hyp = 18 / math.sqrt(5)

# For calculation clarity, we take the one that gives higher value
dx_hyp_final = max(abs(dx_hyp), abs(dy_hyp))
dy_hyp_final = min(abs(dx_hyp), abs(dy_hyp))

optimal_sides = [
    (dx_leg1, dy_leg1),      # Leg 1
    (dy_leg2, dx_leg2),      # Leg 2, perpendicular, so projections swap
    (dx_hyp_final, dy_hyp_final)  # Hypotenuse
]

# We need absolute values for Δx, Δy
abs_optimal_sides = [(abs(s[0]), abs(s[1])) for s in optimal_sides]


largest_k = calculate_k(abs_optimal_sides)
print("\nFinal Answer:")
print("The largest number k is:", largest_k)
