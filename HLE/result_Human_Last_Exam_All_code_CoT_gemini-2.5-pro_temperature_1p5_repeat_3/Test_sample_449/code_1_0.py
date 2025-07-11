import math

def solve():
    """
    Calculates the approximate probability that the conditioned random walk
    never enters the set of the four neighbours of the origin.
    """
    # The starting point is (3000, 4000)
    x, y = 3000, 4000

    # Calculate the distance from the origin
    distance = math.sqrt(x**2 + y**2)

    # The probability of hitting the set of neighbors is approximately 2 / ln(distance)
    # The derivation is explained in the text above.
    prob_hit = 2 / math.log(distance)

    # The probability of never hitting the set is 1 - prob_hit
    prob_never_hit = 1 - prob_hit

    # Output the equation and the result
    print(f"The probability is calculated by the equation: 1 - 2 / ln({int(distance)})")
    print(f"1 - 2 / {math.log(distance):.4f} = {prob_never_hit:.4f}")
    
    # Give the approximate answer with two significant digits
    approx_answer = f"{prob_never_hit:.2g}"
    print(f"The approximate answer with two significant digits is: {approx_answer}")
    
    # Final answer format
    print(f"<<<{approx_answer}>>>")

solve()
