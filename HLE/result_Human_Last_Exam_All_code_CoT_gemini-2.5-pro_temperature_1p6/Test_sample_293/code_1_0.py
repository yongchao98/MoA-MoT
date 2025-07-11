import math

# The value of the integral is given by the formula pi * (1 - log(2*pi))
pi_val = math.pi
integral_val = pi_val * (1 - math.log(2 * pi_val))

# The integral was computed as I = Integral(log(sin(x))dx) - Integral(log(x)dx) from 0 to pi.
# The final equation is I = (-pi*log(2)) - (pi*log(pi) - pi) = pi - pi*log(2) - pi*log(pi).
# Let's print the values of the terms in the final expression.
term_pi = pi_val
term_pi_log2 = pi_val * math.log(2)
term_pi_logpi = pi_val * math.log(math.pi)
final_equation_str = f"pi - pi*log(2) - pi*log(pi)"
computed_value = term_pi - term_pi_log2 - term_pi_logpi

print(f"The integral evaluates to the expression: {final_equation_str}")
print(f"The numerical components are:")
print(f"pi = {term_pi}")
print(f"pi*log(2) = {term_pi_log2}")
print(f"pi*log(pi) = {term_pi_logpi}")
print(f"Final computation: {term_pi} - {term_pi_log2} - {term_pi_logpi} = {computed_value}")
print(f"The final value of the integral is: {integral_val}")
