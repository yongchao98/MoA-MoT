import math

def solve():
    """
    This function calculates the maximum number of grid squares the triangle's
    perimeter can cross.
    
    The number of squares 'k' can be expressed as a function of the triangle's 
    orientation angle, theta. For 0 <= theta <= 45 degrees, the formula is:
    k(theta) = 4 * ceil(18*cos(theta)) + 2 * ceil(18*sin(theta))
    
    We find the maximum value of this function by testing critical points.
    """
    
    max_k = 0
    best_u = 0
    
    # We test values for u = 18*sin(theta) that are just above an integer I.
    # u ranges from 0 up to 9*sqrt(2) which is approx 12.72.
    # So we test for integers I from 0 to 12.
    
    # Test when u = I + epsilon
    for i in range(13):
        u = i + 0.000001
        if u * u >= 324:
            continue
        
        cos_val = math.sqrt(324 - u * u)
        current_k = 4 * math.ceil(cos_val) + 2 * math.ceil(u)
        
        if current_k > max_k:
            max_k = current_k
            best_u = u

    # For completeness, also check when v = 18*cos(theta) is just above an integer M.
    # v ranges from 18 down to ~12.72. So we test for M from 18 down to 13.
    for m in range(18, 12, -1):
        v = m + 0.000001
        if v * v >= 324:
            continue

        u = math.sqrt(324 - v * v)
        # We must be in the theta <= 45 deg range, which means u <= v
        if u > v:
            continue
        
        current_k = 4 * math.ceil(v) + 2 * math.ceil(u)

        if current_k > max_k:
            max_k = current_k
            # No need to update best_u, as it's just for reporting.
            
    # Now, output the final equation details for the found maximum.
    u = best_u
    cos_val_at_max = math.sqrt(324 - u*u)
    
    ceil_cos = math.ceil(cos_val_at_max)
    ceil_sin = math.ceil(u)
    
    final_k = 4 * ceil_cos + 2 * ceil_sin

    print("The maximum value for k is found with a specific orientation.")
    print("The components of the calculation for the maximum k are:")
    print(f"  ceil(18 * cos(theta)) = {ceil_cos}")
    print(f"  ceil(18 * sin(theta)) = {ceil_sin}")
    print("\nThe final equation is:")
    print(f"k = 4 * {ceil_cos} + 2 * {ceil_sin}")
    print(f"k = {4 * ceil_cos} + {2 * ceil_sin}")
    print(f"k = {final_k}")

solve()
<<<86>>>