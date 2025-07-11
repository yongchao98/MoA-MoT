import math

def solve():
    """
    Calculates the expected number of hits on yellow circles based on geometric analysis.
    """
    # Step 1: Define given measurements and find the radius of the large white circles.
    r_yellow = 1  # cm
    
    # We found that for integer coordinates, the smallest radius for the large
    # white circles is R=4.
    R_white = 4 # cm
    
    # Step 2: Determine the target's height H and width W based on integer constraints.
    # From the Pythagorean triple (6,8,10) for staggering and middle row rule:
    H = 24 # cm
    
    # Assuming the width of a green bar w_g is equal to R_white:
    w_g = R_white
    W = 2 * w_g + 7 * R_white
    
    # Step 3: Calculate the total area of the target.
    total_area = W * H
    
    # Step 4: Calculate the total area of the yellow circles.
    num_yellow_circles = 6
    area_yellow_circle = math.pi * (r_yellow ** 2)
    total_yellow_area = num_yellow_circles * area_yellow_circle
    
    # Step 5: Calculate the expected number of hits out of 10000 shots.
    num_shots = 10000
    expected_hits = (total_yellow_area / total_area) * num_shots
    
    # Print the equation
    print("Equation for expected hits:")
    print(f"Expected Hits = (Number of Yellow Circles * Area of one Yellow Circle) / (Target Width * Target Height) * Number of Shots")
    print(f"Expected Hits = ({num_yellow_circles} * pi * {r_yellow}^2) / ({W} * {H}) * {num_shots}")
    print(f"Expected Hits = ({num_yellow_circles} * pi) / {total_area} * {num_shots}")
    print(f"Expected Hits = ({num_shots * num_yellow_circles} * pi) / {total_area}")
    print(f"Expected Hits = {num_shots * num_yellow_circles / total_area:.4f} * pi")

    print(f"\nFinal calculation:")
    print(f"Expected Hits = {expected_hits}")
    
    # The final answer needs to be returned in a specific format
    # The calculated value is ~218.17
    # Depending on rounding rules (e.g. to nearest integer), it could be 218.
    # The question doesn't specify, so let's provide the unrounded value.
    return expected_hits

result = solve()
# We will present the answer rounded to one decimal place as it is a statistical expectation.
# final_answer = round(result, 1)
# Or rounding to the nearest integer
final_answer = round(result)
print(f"\nRounded to the nearest integer, the expected number of hits is {final_answer}.")
# For the required format, let's output the value directly
# The result is ~218.166... so rounding to one decimal place is 218.2

# The format wants a single number. 
# print(f'<<<{round(result, 1)}>>>')
# Let's provide the value rounded to the nearest integer as is common for counting "times".
# print(f'<<<{final_answer}>>>')

# I'll stick to one decimal place since it's an expectation value.
# print(f"<<<{round(result,1)}>>>")
# No, let's output the result from the variable directly as requested in other similar problems
# I will use my final integer rounded answer
final_answer_to_output = 218