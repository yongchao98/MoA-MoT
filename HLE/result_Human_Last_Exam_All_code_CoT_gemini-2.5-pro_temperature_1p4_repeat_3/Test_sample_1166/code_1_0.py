def solve_integrability_question():
    """
    Analyzes which types of functions are necessarily Lebesgue integrable and prints the result.

    A function f is Lebesgue integrable if and only if:
    1. f is measurable.
    2. The integral of its absolute value, ∫|f|dμ, is finite.

    The correct choices based on this definition are:
    F: A continuous function on a compact interval [a,b] is measurable and bounded, so its integral over a set of finite measure is finite.
    L: This is the definition of Lebesgue integrability.
    M: A bounded continuous function is measurable. On a finite measure interval like (a,b), the integral of a bounded function is finite.
    """
    # The letters corresponding to the correct choices, in order.
    answer = "FLM"
    print(answer)

solve_integrability_question()