import math

def solve():
    """
    Solves for the Lebesgue measure of the set S multiplied by 10^6.
    
    The set S consists of all starting points x_0 in [0,1] for which the sequence 
    x_{n+1} = f(x_n) with f(x) = (2x + sin(2*pi*x))/3 contains exactly 7 distinct values.

    1. A sequence with a finite number of distinct values must be eventually periodic.
    2. This means that for any x_0 in S, x_0 is an eventually periodic point of the map f.
    3. The set of eventually periodic points for a C^1 map on an interval (that is not piecewise constant)
       is a countable set. This is because the set of periodic points of any period p is finite,
       and the set of all preimages of a countable set is also countable.
    4. A countable set has a Lebesgue measure of zero.
    5. The set S is a subset of the set of all eventually periodic points, so its Lebesgue measure
       is also zero.
    6. Therefore, the Lebesgue measure of S multiplied by 10^6 is 0.
    """
    
    # Based on the mathematical reasoning, the Lebesgue measure of S is 0.
    lebesgue_measure_of_S = 0
    
    # The problem asks for this measure multiplied by 10^6.
    result = lebesgue_measure_of_S * 10**6
    
    print(int(result))

solve()