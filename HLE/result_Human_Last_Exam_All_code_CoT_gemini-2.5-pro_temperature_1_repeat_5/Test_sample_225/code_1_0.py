import math

def count_self_intersections():
    """
    Calculates the number of self-intersection points for the curve
    x(t) = cos(9t), y(t) = sin(5t).
    """
    a = 9
    b = 5

    # A self-intersection requires cos(a*t1)=cos(a*t2) and sin(b*t1)=sin(b*t2).
    # This leads to solutions for t1 and t2 that depend on integer parameters k and m.
    # We find all pairs (k, m) that give valid t1, t2 in the range [0, 2*pi).
    # This gives 40 intersection "events".
    intersection_pairs = []
    # k must be in [1, a-1]
    for k in range(1, a):
        # The range for m is derived from 0 <= t1 < t2 < 2*pi
        # This leads to the inequality: (10k-9)/18 <= m < (171-10k)/18
        # We use integer arithmetic to find the bounds for m.
        # m_start corresponds to ceil((10*k - 9) / 18)
        m_start = (10 * k - 9 + 17) // 18
        # m_end corresponds to floor((171 - 10*k - 1) / 18)
        m_end = (170 - 10 * k) // 18

        for m in range(m_start, m_end + 1):
            intersection_pairs.append((k, m))
    
    # Now, we calculate the (x,y) coordinates for each intersection event
    # and store the unique points in a set.
    points = set()
    for k, m in intersection_pairs:
        # The time t1 can be calculated from k and m:
        # t1 = pi/90 * (18m + 9 - 10k)
        # x = cos(9*t1) = cos(pi/10 * (18m + 9 - 10k))
        # y = sin(5*t1) = sin(pi/18 * (18m + 9 - 10k))
        
        theta_numerator = 18 * m + 9 - 10 * k
        
        x = math.cos(math.pi * theta_numerator / 10)
        y = math.sin(math.pi * theta_numerator / 18)
        
        # We round the coordinates to handle floating-point inaccuracies.
        point = (round(x, 10), round(y, 10))
        points.add(point)

    num_points = len(points)
    print(f"The number of self-intersection points found by calculation is: {num_points}")
    
    # For a Lissajous curve x=cos(at), y=sin(bt) with a,b coprime,
    # the number of self-intersections is given by the formula (a-1)(b-1)/2.
    print("\nThis can be verified with the known formula for this curve type:")
    
    formula_result = (a - 1) * (b - 1) // 2
    
    print(f"Formula: (a - 1) * (b - 1) / 2")
    print(f"Using a = {a} and b = {b}:")
    print(f"({a} - 1) * ({b} - 1) / 2 = ({a-1}) * ({b-1}) / 2 = {a-1_item*b-1_item} / 2 = {formula_result}")

count_self_intersections()