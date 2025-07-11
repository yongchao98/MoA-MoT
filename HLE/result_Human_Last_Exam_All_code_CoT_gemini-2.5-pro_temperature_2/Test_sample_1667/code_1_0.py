import numpy as np

def in_circle_test(p1, p2, p3, p4):
    """
    Tests if point p4 is inside the circumcircle of points p1, p2, p3.
    This uses a determinant test. The sign of the determinant indicates the position.
    A positive result here means p4 is inside the circle defined by p1,p2,p3,
    assuming p1,p2,p3 are ordered counter-clockwise.
    If the order is random, the absolute value is the same, but the sign might flip.
    However, for the purpose of this simulation, we check the sign of two
    determinants. If D(p1,p2,p3,p4) and D(p1,p2,p4,p3) have opposite signs,
    then the logic holds, and we have a 50% chance for a specific point.
    For simplicity, we compute one determinant. A positive result corresponds
    to one of the two pairs from the "2-in, 2-out" property.

    Args:
        p1, p2, p3, p4: Points as numpy arrays, e.g., np.array([x, y]).
    Returns:
        A float representing the determinant. >0 if inside (for CCW points), <0 if outside.
    """
    # The matrix for the in-circle test.
    # It checks the orientation of the lifted points in 3D.
    matrix = np.array([
        [p1[0], p1[1], p1[0]**2 + p1[1]**2, 1],
        [p2[0], p2[1], p2[0]**2 + p2[1]**2, 1],
        [p3[0], p3[1], p3[0]**2 + p3[1]**2, 1],
        [p4[0], p4[1], p4[0]**2 + p4[1]**2, 1]
    ])

    # To ensure the test works regardless of orientation (clockwise vs counter-clockwise)
    # we can check against the orientation of the three base points.
    orientation_matrix = np.array([
        [p1[0], p1[1], 1],
        [p2[0], p2[1], 1],
        [p3[0], p3[1], 1]
    ])
    
    det = np.linalg.det(matrix)
    orientation_det = np.linalg.det(orientation_matrix)

    # We return the determinant value, but normalized by the triangle orientation.
    # This ensures a positive result means "inside" consistently.
    return det * orientation_det

def run_simulation(num_trials):
    """
    Runs a Monte Carlo simulation to estimate the probability.
    """
    inside_circle_count = 0

    for _ in range(num_trials):
        # Generate 4 random points (ducks) in a 1x1 unit square
        points = np.random.rand(4, 2)
        p1, p2, p3, p4 = points[0], points[1], points[2], points[3]
        
        # We need to handle the (near-zero probability) case of p1,p2,p3 being collinear.
        # In this case, the circumcircle is undefined (a line). We can just skip such trials.
        # The orientation determinant will be close to zero.
        orientation_det_check = np.linalg.det(np.array([[p1[0], p1[1], 1], [p2[0], p2[1], 1], [p3[0], p3[1], 1]]))
        if abs(orientation_det_check) < 1e-9:
            continue

        # Check if the 4th duck is inside the circle of the first three
        # A positive result from our test indicates the point is inside.
        if in_circle_test(p1, p2, p3, p4) > 0:
            inside_circle_count += 1
            
    # The probability is the number of successful events divided by the total number of trials
    probability = inside_circle_count / num_trials
    
    print(f"Number of trials: {num_trials}")
    print(f"Times 4th duck was in the circle: {inside_circle_count}")
    print(f"Final calculated probability: {probability}")
    

if __name__ == '__main__':
    # Set the number of simulations to run
    # A higher number will give a more accurate result
    trials = 100000
    run_simulation(trials)
