def solve():
    """
    This function determines which of the given function descriptions are 
    necessarily Lebesgue integrable.
    """
    
    # A list of letters corresponding to the properties that guarantee
    # a function is Lebesgue integrable, based on the analysis.
    # F: A continuous function on [a,b] is measurable and bounded, so its integral is finite.
    # L: This is the definition of Lebesgue integrability (measurable and |f| has a finite integral).
    # M: A bounded continuous function on (a,b) is measurable, and its integral over a finite-measure
    #    domain is finite.
    integrable_choices = ['F', 'L', 'M']
    
    # The final answer is a string of these letters in order.
    answer = "".join(integrable_choices)
    
    print(answer)

solve()