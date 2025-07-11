import math

def calculate_cube_curve_percentage():
    """
    Solves the cube curve problem by calculating the total length of the locus C
    and expressing its ratio to 2*pi*r as a percentage.
    """
    print("This problem asks to calculate a ratio related to the length of a curve C on a cube's surface.")
    print("The final result is independent of the cube's side length r.")

    # The curve C is found to be composed of several distinct circular arcs on the cube's faces.
    # By unfolding the cube, we can determine the number and length of these arcs.
    # It turns out all elementary arcs have the same length.

    # 1. Determine the number of arcs.
    # - 2 arcs are on the two faces containing the midpoint P.
    # - 4 arcs are on the two faces adjacent to both of the primary faces (2 arcs on each).
    arc_count = 6
    print(f"\nThe curve C is composed of a total of {arc_count} elementary arcs.")

    # 2. Determine the length of each arc.
    # The length of each arc, determined from the geometry on the unfolded cube, is (r * pi / 3).
    arc_len_numerator = 1
    arc_len_denominator = 3
    print(f"The length of each arc is equal to ({arc_len_numerator}/{arc_len_denominator})*pi*r.")

    # 3. Calculate the total length of C by summing the arc lengths.
    total_len_numerator_factor = arc_count * arc_len_numerator
    total_len_denominator_factor = arc_len_denominator
    # This simplifies to (6/3)*pi*r = 2*pi*r
    final_len_factor = int(total_len_numerator_factor / total_len_denominator_factor)

    print(f"\nThe total length of curve C is the number of arcs multiplied by the length of each arc:")
    print(f"Length C = {arc_count} * ({arc_len_numerator}/{arc_len_denominator})*pi*r = {final_len_factor}*pi*r.")

    # 4. Calculate the final ratio and percentage.
    # The problem requires dividing the length of C by 2*pi*r.
    denominator_in_ratio = 2
    ratio = final_len_factor / denominator_in_ratio
    percentage = int(ratio * 100)

    print(f"\nFinally, we divide this total length by {denominator_in_ratio}*pi*r and express it as a percentage.")
    print(f"The calculation is: (({final_len_factor}*pi*r) / ({denominator_in_ratio}*pi*r)) * 100")
    print(f"This simplifies to: ({final_len_factor}/{denominator_in_ratio}) * 100 = {int(ratio)} * 100 = {percentage}%")


calculate_cube_curve_percentage()
print("\n<<<100>>>")