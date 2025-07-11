import math
import random

def get_area_in_square(nx, ny, num_samples_per_square):
    """
    Calculates the area of the region R within a single unit square
    [nx, nx+1) x [ny, ny+1) using Monte Carlo integration.
    """
    # For a point in this square, floor(x)=nx and floor(y)=ny
    k_squared = nx*nx + ny*ny
    
    # We are only calling this function for pairs where k is an integer
    k = int(math.sqrt(k_squared))
    
    # The condition for R simplifies to k <= |z| < k+1
    # Squaring gives k^2 <= x^2+y^2 < (k+1)^2
    count = 0
    inner_radius_sq = k * k
    outer_radius_sq = (k + 1) * (k + 1)

    for _ in range(num_samples_per_square):
        x = random.uniform(nx, nx + 1)
        y = random.uniform(ny, ny + 1)
        r_sq = x*x + y*y
        if inner_radius_sq <= r_sq < outer_radius_sq:
            count += 1

    # The area of the unit square is 1.
    area = float(count) / num_samples_per_square
    return area

def solve():
    """
    Solves the problem by summing up the area contributions from all relevant squares.
    """
    # Pairs (nx, ny) where 0<=nx,ny<=5 and sqrt(nx^2+ny^2) is an integer.
    # We group symmetric pairs like (1,0) and (0,1) together for clarity.
    pairs = [
        (0, 0),
        (1, 0), (0, 1),
        (2, 0), (0, 2),
        (3, 0), (0, 3),
        (4, 0), (0, 4),
        (5, 0), (0, 5),
        (3, 4),
        (4, 3)
    ]

    num_samples_per_square = 3000000  # Number of random points per unit square
    
    area_contributions = []
    
    for nx, ny in pairs:
        area_contribution = get_area_in_square(nx, ny, num_samples_per_square)
        area_contributions.append(area_contribution)

    total_area = sum(area_contributions)

    # Print the full equation showing the sum of each contribution
    equation_parts = [f"{area:.4f}" for area in area_contributions]
    equation_str = " + ".join(equation_parts)
    
    print(f"The area of R is the sum of contributions from each valid unit square.")
    print("The contributing areas are calculated for squares (nx,ny):")
    print("(0,0), (1,0), (0,1), (2,0), (0,2), (3,0), (0,3), (4,0), (0,4), (5,0), (0,5), (3,4), (4,3)")
    print("\nThe final sum is:")
    print(f"{equation_str} = {total_area:.4f}")
    
    print(f"\nFinal Area rounded to two decimals: {total_area:.2f}")


if __name__ == '__main__':
    solve()