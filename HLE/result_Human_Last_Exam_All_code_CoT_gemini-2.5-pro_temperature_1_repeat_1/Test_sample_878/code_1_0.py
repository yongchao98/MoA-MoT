import math

def solve_egg_drop():
    """
    Calculates the egg drop time using simple integer and fractional arithmetic.
    """
    print("Yes, your son can calculate the time. Here is a simple way:")
    print("-" * 50)

    # Step 1: Find the height of the skyscraper (h)
    print("Step 1: Find the height (h).")
    print("The angle is 1/4 of a 90 degree angle, and the distance is 240m.")
    print("We can approximate tan(22.5 degrees) with the simple fraction 5/12.")

    d = 240
    tan_alpha_num = 5
    tan_alpha_den = 12
    h = int(d * tan_alpha_num / tan_alpha_den)

    print(f"h ≈ {d} * ({tan_alpha_num} / {tan_alpha_den}) = {h} m")
    print("-" * 50)

    # Step 2: Find the time of fall (t)
    print("Step 2: Find the time to fall (t).")
    print("The formula is t = sqrt(2 * h / g). We will use g ≈ 10 m/s^2.")
    print("Here is the calculation with all the numbers:")

    g = 10
    two = 2
    val_inside_sqrt_1 = int(two * h)
    val_inside_sqrt_2 = int((two * h) / g)
    four = 4
    five = 5
    sqrt5_approx_num = 9
    sqrt5_approx_den = 4
    final_num = two * sqrt5_approx_num
    final_den = sqrt5_approx_den
    final_result = final_num / final_den

    print(f"t ≈ sqrt( ({two} * {h}) / {g} )")
    print(f"  = sqrt( {val_inside_sqrt_1} / {g} )")
    print(f"  = sqrt( {val_inside_sqrt_2} )")
    print(f"  = sqrt( {four} * {five} )")
    print(f"  = {two} * sqrt({five})")
    print(f"  ≈ {two} * ({sqrt5_approx_num}/{sqrt5_approx_den})")
    print(f"  = {final_num}/{final_den}")
    print(f"  = {final_result} seconds")
    print("-" * 50)

solve_egg_drop()
# The largest integer used in the calculation is 240.
# The calculation is possible with an error < 0.1s.
# So the answer is Yz where z=240.
print("<<<Y240>>>")
