# This script analyzes the properties of functions to determine if they
# are necessarily Lebesgue integrable.

# A function f is Lebesgue integrable if it meets two conditions:
# 1. f must be a measurable function.
# 2. The integral of its absolute value must be finite (∫|f|dμ < ∞).

# We evaluate each choice against these two conditions.

# A. A bounded function: Fails (1). Not always measurable.
# B. A bounded measurable function: Fails (2) on ℝ. e.g., f(x) = 1.
# C. A measurable function: Fails (2) on ℝ. e.g., f(x) = x.
# D. A continuous function: Fails (2) on ℝ. e.g., f(x) = x.
# E. A measurable function on [a,b]: Fails (2). e.g., f(x) = 1/x on [0,1].
# F. A continuous function on [a,b]: Correct. Is measurable and bounded on [a,b], so ∫|f| is finite.
# G. A bounded function on [a,b]: Fails (1). Not always measurable.
# H. A function whose absolute value is integrable: Fails (1). f itself might not be measurable.
# I. A function whose absolute value is integrable on [a,b]: Fails (1), same as H.
# J. A continuous function on (a,b): Fails (2). e.g., f(x) = 1/(x-a).
# H. (duplicate) A bounded function on (a,b): Fails (1). Not always measurable.
# K. A measurable function on (a,b): Fails (2). e.g., f(x) = 1/(x-a).
# L. A measurable function whose absolute value is integrable: Correct. This is the definition of Lebesgue integrability.
# M. A bounded continuous function on (a,b): Correct. Is measurable and bounded on a finite interval, so ∫|f| is finite.

# Concatenating the letters of the correct choices in order gives the answer.
correct_choices = "FLM"

print(correct_choices)