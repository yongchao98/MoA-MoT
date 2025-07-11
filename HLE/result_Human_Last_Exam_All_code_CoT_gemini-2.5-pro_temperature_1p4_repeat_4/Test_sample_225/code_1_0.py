def solve_lissajous_intersections():
    """
    Calculates the number of self-intersection points for the Lissajous curve
    x(t) = cos(at), y(t) = sin(bt) where a and b are odd coprime integers.
    """
    # Parameters from the curve equations x(t) = cos(9t), y(t) = sin(5t)
    a = 9
    b = 5

    print(f"The curve is a Lissajous curve defined by x(t) = cos({a}t) and y(t) = sin({b}t).")
    print("For this type of curve, with 'a' and 'b' being coprime odd integers,")
    print("the number of self-intersection points can be found with the formula: (a - 1) * (b - 1) / 2.")
    print("-" * 20)
    print(f"Given parameters: a = {a}, b = {b}")

    # Step-by-step calculation
    a_minus_1 = a - 1
    b_minus_1 = b - 1
    product = a_minus_1 * b_minus_1
    result = product // 2  # Use integer division

    # Print the full equation as requested
    print("The calculation is:")
    print(f"({a} - 1) * ({b} - 1) / 2 = {a_minus_1} * {b_minus_1} / 2 = {product} / 2 = {result}")
    print("-" * 20)
    print(f"The number of self-intersection points is: {result}")

if __name__ == "__main__":
    solve_lissajous_intersections()
