def solve_lissajous_intersections():
    """
    Calculates the number of self-intersection points for the curve
    (cos(at), sin(bt)) for given a and b.
    """
    # The parameters for the curve x(t) = cos(9t), y(t) = sin(5t)
    a = 9
    b = 5

    print(f"The curve is a Lissajous curve with parameters a = {a} and b = {b}.")
    print("To find the number of self-intersections, we can use a standard formula.")
    
    # The formula applies if a and b are coprime and both are odd.
    # Let's verify the conditions.
    # Condition 1: a and b are coprime (GCD is 1).
    p, q = a, b
    while q:
        p, q = q, p % q
    gcd = p
    
    # Condition 2: a and b are both odd.
    are_both_odd = (a % 2 != 0) and (b % 2 != 0)

    if gcd == 1 and are_both_odd:
        print(f"The parameters a={a} and b={b} are coprime and both are odd.")
        print("The formula for the number of self-intersections is: N = (a - 1) * (b - 1) / 2")
        
        # Perform the calculation step-by-step
        a_minus_1 = a - 1
        b_minus_1 = b - 1
        numerator = a_minus_1 * b_minus_1
        num_intersections = numerator / 2
        
        print("\nCalculating the result:")
        print(f"N = ({a} - 1) * ({b} - 1) / 2")
        print(f"N = {a_minus_1} * {b_minus_1} / 2")
        print(f"N = {numerator} / 2")
        print(f"N = {int(num_intersections)}")
        
    else:
        print("The conditions for the simple formula are not met.")

solve_lissajous_intersections()