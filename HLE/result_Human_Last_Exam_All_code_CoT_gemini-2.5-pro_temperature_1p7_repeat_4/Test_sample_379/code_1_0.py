import math

def solve_cube_locus_problem():
    """
    Solves the cube locus problem as described.
    
    This function calculates the length of a locus of points C on a cube's surface,
    divides it by 2*pi*r, and expresses the result as a whole number percentage.
    """
    
    # We can assume r=1, as the value 'r' will cancel out in the final ratio.
    r = 1
    pi = math.pi
    
    print("Step 1: Define constants (we can use r=1 as it's a ratio).")
    print(f"r = {r}")
    print(f"pi = {pi:.4f}\n")
    
    # Step 2: Calculate the length of one arc.
    # The locus consists of 4 identical arcs on the 4 faces adjacent to P's edge.
    # The length of each arc is (pi * r) / 3.
    num_arcs = 4
    length_one_arc = (pi * r) / 3
    
    print(f"Step 2: The locus C consists of {num_arcs} identical arcs.")
    print(f"The length of a single arc is (pi * r) / 3 = {length_one_arc:.4f}\n")

    # Step 3: Calculate the total length of the locus C.
    total_length_C = num_arcs * length_one_arc

    print("Step 3: Calculate the total length of C.")
    print(f"Total Length = {num_arcs} * {length_one_arc:.4f} = {total_length_C:.4f}\n")
    
    # Step 4: Calculate the divisor 2 * pi * r.
    divisor = 2 * pi * r

    print("Step 4: Calculate the divisor 2 * pi * r.")
    print(f"Divisor = 2 * {pi:.4f} * {r} = {divisor:.4f}\n")
    
    # Step 5: Compute the ratio.
    ratio = total_length_C / divisor
    
    print("Step 5: Compute the ratio of the total length to the divisor.")
    print(f"Ratio = {total_length_C:.4f} / {divisor:.4f} = {ratio:.4f}\n")
    
    # Step 6: Convert to a whole number percentage.
    percentage = round(ratio * 100)
    
    print("Step 6: Convert the ratio to a whole number percentage and show the final equation.")
    # We display the numbers in the final calculation as requested.
    print(f"Final Answer (%) = round( (Total Length / Divisor) * 100 )")
    print(f"                 = round( ( {total_length_C:.3f} / {divisor:.3f} ) * 100 )")
    print(f"                 = round( {ratio:.3f} * 100 )")
    print(f"                 = {percentage}\n")

solve_cube_locus_problem()

# The final answer in the required format
final_answer = round(( (4/3) * math.pi / (2 * math.pi) ) * 100)
print(f"<<<{final_answer}>>>")