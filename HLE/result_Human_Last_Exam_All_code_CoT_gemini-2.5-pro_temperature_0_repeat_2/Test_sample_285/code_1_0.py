import sys

def solve_integrability():
    """
    Calculates the largest p such that the function I is not in L^p(R^9).
    """
    # The dimension of the parameter space a = (a1, ..., a9) is d.
    d = 9

    # The slowest decay rate of the integral I(a) behaves as |a|^(-delta)
    # as |a| -> infinity.
    # This slowest decay occurs when the phase polynomial is made as "flat" as possible.
    # By choosing the parameters a_i, the phase polynomial can be made to be of the form,
    # for example, lambda * x^3.
    # The integral becomes integral(exp(2*pi*i * lambda * x^3)) dx dy,
    # which decays as lambda^(-1/3).
    # Therefore, the minimal decay exponent is delta = 1/3.
    delta_numerator = 1
    delta_denominator = 3
    
    # The function I(a) is in L^p(R^9) if the integral of |I(a)|^p converges.
    # This happens if p * delta > d.
    # It is not in L^p(R^9) if p * delta <= d.
    # We want the largest p such that I is not in L^p.
    # This threshold is given by p = d / delta.
    p = d * delta_denominator / delta_numerator

    print(f"The dimension of the parameter space is d = {d}.")
    print(f"The minimal decay exponent of the integral is delta = {delta_numerator}/{delta_denominator}.")
    print(f"The largest p such that the function I is not in L^p(R^9) is given by the equation p = d / delta.")
    print(f"p = {d} / ({delta_numerator}/{delta_denominator}) = {int(p)}")

solve_integrability()