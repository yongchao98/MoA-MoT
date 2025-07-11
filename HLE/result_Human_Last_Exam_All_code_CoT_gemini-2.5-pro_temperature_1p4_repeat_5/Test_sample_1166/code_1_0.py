# This script prints the letters corresponding to the choices that are 
# necessarily Lebesgue integrable.
# The analysis is based on the definition of Lebesgue integrability, which requires a function
# to be both measurable and have a finite integral for its absolute value.

# F: A continuous function on a closed interval [a,b] is bounded, and a bounded
#    measurable function on a finite measure space is integrable.
# L: This is the definition of a Lebesgue integrable function.
# M: A bounded continuous function on (a,b) is a bounded measurable function
#    on a space of finite measure, which is integrable.

correct_choices = "FLM"
print(correct_choices)