import math

def solve_chain_problem():
    """
    This function calculates the time when the chain first loses contact with the ground.
    """
    # Given parameters
    h = 1.0  # m
    r_cm = 25.0  # cm
    r = r_cm / 100.0  # m
    l_c = 10.0  # m
    v = 10.0  # m/s
    d = 20.0  # m
    l_shadow = 10.0 * math.sqrt(3)  # m
    beta_dot = 1.0  # rad/s

    # Step 1: Calculate geometric parameters
    R = d / 2.0
    # From l = d * cos(alpha)
    cos_alpha = l_shadow / d
    alpha = math.acos(cos_alpha)
    sin_alpha = math.sin(alpha)

    # Robot angular speed along the path
    omega = v / R

    # Since omega = 1 rad/s and beta_dot = 1 rad/s, the time 't' is equal
    # to the angle of robot travel and the angle of arm rotation.

    # Step 2: Set up the quadratic equation for s = sin(t)
    # The equation is of the form: A*s^2 + B*s + C = 0
    # A = -r*sin(alpha)
    # B = R*sin(alpha) + r*cos(alpha)
    # C = R*sin(alpha) + h*cos(alpha) + r*sin(alpha) - l_c

    A = -r * sin_alpha
    B = R * sin_alpha + r * cos_alpha
    C = (R * sin_alpha) + (h * cos_alpha) + (r * sin_alpha) - l_c

    print(f"The final equation to solve for s=sin(t) is of the form A*s^2 + B*s + C = 0.")
    print(f"The coefficients are:")
    print(f"A = -r*sin(alpha) = {A}")
    print(f"B = R*sin(alpha) + r*cos(alpha) = {B}")
    print(f"C = R*sin(alpha) + h*cos(alpha) + r*sin(alpha) - l_c = {C}")
    print("\nSo the final equation is:")
    print(f"({A})*sin(t)^2 + ({B})*sin(t) + ({C}) = 0")


    # Step 3: Solve the quadratic equation
    discriminant = B**2 - 4*A*C

    if discriminant < 0:
        print("\nNo real solution for sin(t) exists. The chain never leaves the ground under these conditions.")
        s = None
    else:
        # We find the two possible solutions for s = sin(t)
        s1 = (-B + math.sqrt(discriminant)) / (2 * A)
        s2 = (-B - math.sqrt(discriminant)) / (2 * A)

        # Check which solution is valid (-1 <= s <= 1)
        valid_s = []
        if -1 <= s1 <= 1:
            valid_s.append(s1)
        if -1 <= s2 <= 1:
            valid_s.append(s2)

        if not valid_s:
            print(f"\nSolutions for sin(t) are s1={s1}, s2={s2}. Both are outside the valid range [-1, 1].")
            s = None
        else:
            # We want the first time, so we look for the smallest positive t.
            # This will come from the smallest positive s.
            s = min(val for val in valid_s if val >= 0) if any(val >= 0 for val in valid_s) else None
            if s is None:
                # If both valid_s are negative, this means t is > pi. Take smallest |s|.
                s = max(valid_s)
    
    if s is not None:
      # Calculate the time t
      t = math.asin(s)
      # We are looking for the first positive time. asin returns in [-pi/2, pi/2].
      # If t is negative, the first positive time is in the next cycle, but the
      # problem implies a solution is found quickly.
      # A positive sin(t) gives a t in [0, pi/2], which is the earliest possible positive time.
      if t < 0:
          # For a negative sin(t), asin gives t in [-pi/2, 0]. The first positive time would be pi - t
          t = math.pi - t

      print(f"\nThe valid solution for sin(t) is: {s}")
      print(f"\nThe time when the chain first loses contact with the ground is t = arcsin({s:.4f}) seconds.")
      print(f"t = {t}")
      return t
    else:
      return None

result = solve_chain_problem()