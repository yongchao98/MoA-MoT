def solve():
    """
    Analyzes the properties of functions to determine if they are necessarily Lebesgue integrable.
    """
    # Based on the analysis:
    # F. A continuous function on [a,b] is integrable.
    # L. A measurable function whose absolute value is integrable is the definition of integrable.
    # M. A bounded continuous function on (a,b) is integrable.
    # The other options are not necessarily integrable due to potential non-measurability,
    # unboundedness, or an infinite measure of the domain.
    answer = "FLM"
    print(answer)

solve()