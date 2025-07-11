import math
from fractions import Fraction

def solve_egg_drop_problem():
    """
    This function calculates the egg drop time using simple integer and
    fraction arithmetic, checks if the accuracy is within the required
    bounds, and identifies the largest integer used in the calculation.
    """
    # --- GIVEN VALUES AND APPROXIMATIONS ---
    # The problem provides the distance d.
    d = 240

    # The angle is 22.5 degrees. We approximate tan(22.5) as 5/12.
    tan_alpha_num = 5
    tan_alpha_den = 12

    # We use a simple approximation for gravity, g = 10 m/s^2.
    g = 10

    # During the calculation, we'll need to approximate sqrt(5). We'll use 9/4.
    sqrt5_num = 9
    sqrt5_den = 4

    # This list will track all integers used in the calculation.
    used_integers = set()

    print("Yes, the son can calculate the time. Here is a possible method using simple numbers:")
    print("-" * 70)

    # --- STEP 1: CALCULATE HEIGHT (h) ---
    # To keep numbers small, we divide before multiplying.
    h_intermediate = d / tan_alpha_den
    h = h_intermediate * tan_alpha_num
    used_integers.update([d, tan_alpha_num, tan_alpha_den, int(h_intermediate), int(h)])

    print("1. Calculate the height of the skyscraper (h):")
    # The prompt requires printing the numbers in the final equation.
    print(f"h = {d} * {tan_alpha_num} / {tan_alpha_den} = {int(h)} m")

    # --- STEP 2: CALCULATE TIME SQUARED (t^2) ---
    t_squared_num = 2 * int(h)
    t_squared_val = t_squared_num / g
    used_integers.update([2, int(h), g, t_squared_num, int(t_squared_val)])

    print("\n2. Use the free fall formula t^2 = (2 * h) / g, with g ≈ 10:")
    print(f"t^2 = (2 * {int(h)}) / {g} = {int(t_squared_val)}")

    # --- STEP 3: SOLVE FOR TIME (t) ---
    # t = sqrt(20) = sqrt(4 * 5) = 2 * sqrt(5)
    factor_out = 2  # from sqrt(4)
    sqrt_val = 5   # the remaining part under the square root
    used_integers.update([int(t_squared_val), factor_out, sqrt_val])
    
    print("\n3. Solve for t:")
    print(f"t = sqrt({int(t_squared_val)}) which is sqrt(4 * {sqrt_val}) = {factor_out} * sqrt({sqrt_val})")
    
    # Apply the approximation for sqrt(5)
    t_final_num = factor_out * sqrt5_num
    t_final_den = sqrt5_den
    t_final_fraction = Fraction(t_final_num, t_final_den)
    t_final_value = float(t_final_fraction)
    used_integers.update([sqrt5_num, sqrt5_den, t_final_num, t_final_den,
                          t_final_fraction.numerator, t_final_fraction.denominator])

    print(f"Using the approximation sqrt({sqrt_val}) ≈ {sqrt5_num} / {sqrt5_den}:")
    print(f"t = {factor_out} * {sqrt5_num} / {sqrt5_den} = {t_final_fraction.numerator} / {t_final_fraction.denominator} s")
    print(f"\nFinal calculated time: {t_final_value} s")
    print("-" * 70)

    # --- ERROR CHECK (Internal verification) ---
    h_actual = 240 * (math.sqrt(2) - 1)
    g_actual = 9.80665  # Standard gravity
    t_actual = math.sqrt(2 * h_actual / g_actual)
    error = abs(t_final_value - t_actual)
    
    can_solve = error < 0.1

    # --- FINAL ANSWER ---
    if can_solve:
        largest_integer = max(used_integers)
        answer = f"Y{largest_integer}"
    else:
        answer = "N0"
    
    print(f"<<<{answer}>>>")

solve_egg_drop_problem()