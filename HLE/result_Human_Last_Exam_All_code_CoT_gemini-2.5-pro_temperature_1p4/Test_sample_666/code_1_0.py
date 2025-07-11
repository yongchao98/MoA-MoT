import numpy as np

def count_enclosed_poles():
    """
    Counts the number of poles of f(z) enclosed by the curve Gamma.

    The equation of the curve Gamma in (X,Y) = (x+y, x-y) coordinates is:
    P(X,Y) = 3(X^2 - 400)^3 + 1200Y^4 - 20X^3Y^2 = 0
    A pole z = x+iy = a + 2*pi*k*i is inside the curve if P(x+y, x-y) < 0.
    """
    
    # We found the curve is bounded by approximately |X| < 23.96 and |Y| < 10.71.
    # This gives us a search range for a and k.
    # |a + 2*pi*k| < 23.96 and |a - 2*pi*k| < 10.71
    # 2*pi is approx 6.28. If |k| >= 3, |Y| is at least | |a|-18.84 |,
    # combined with |X| bound, it's unlikely to have solutions.
    # We search a safe range for k.
    k_range = range(-5, 6)
    a_range = range(-2024, 2025)
    
    pole_count = 0
    
    pole_counts_by_k = {}

    for k in k_range:
        current_k_count = 0
        for a in a_range:
            x = float(a)
            y = 2 * np.pi * k
            
            X = x + y
            Y = x - y
            
            # Efficiently skip points outside the bounding box.
            if abs(X) >= 23.96 or abs(Y) >= 10.71:
                continue

            # The inside of the curve corresponds to P(X,Y) < 0
            # as P(0,0) = 3*(-400)^3 which is negative.
            val = 3 * (X**2 - 400)**3 + 1200 * Y**4 - 20 * X**3 * Y**2
            
            if val < 0:
                current_k_count += 1
        
        if current_k_count > 0:
            pole_counts_by_k[k] = current_k_count
        pole_count += current_k_count

    print("The final integral is given by 2 * pi * i * N, where N is the number of enclosed poles.")
    print("The number of enclosed poles for each value of k is:")
    for k, count in sorted(pole_counts_by_k.items()):
        print(f"k = {k}: {count} poles")

    print(f"\nThe total number of enclosed poles is N = {pole_count}.")
    
    coefficient = 2 * pole_count
    print(f"\nThe final equation is:")
    print(f"oint_Gamma f(z) dz = 2 * pi * i * {pole_count} = {coefficient} * pi * i")

count_enclosed_poles()