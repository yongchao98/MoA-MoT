def solve():
    """
    This function determines which of the given function types are necessarily Lebesgue integrable.
    F. A continuous function on [a,b]: Yes. It's measurable and bounded on a set of finite measure.
    L. A measurable function whose absolute value is integrable: Yes. This is the definition of Lebesgue integrability.
    M. A bounded continuous function on (a,b): Yes. It's a bounded measurable function on a set of finite measure.
    The final answer is a string of the letters corresponding to the correct choices, in order.
    """
    # The letters corresponding to the necessarily Lebesgue integrable functions are F, L, and M.
    answer = "FLM"
    print(answer)

solve()