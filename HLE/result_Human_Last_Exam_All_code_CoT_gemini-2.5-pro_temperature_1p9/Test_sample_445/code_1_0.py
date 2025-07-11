# The problem is to find the maximal probability p of success for Alice.
# Our analysis concluded that Alice's optimal strategy is to open as many boxes as possible.
# The number of boxes in total is 20. Alice cannot open all of them.
# So, the optimal number of boxes for Alice to open is k=19.

k = 19

# The strategy is to open k=19 boxes and observe the values inside.
# Let the minimum observed value be m and the maximum be M.
# Alice then guesses that the number in the single closed box lies within the interval [m, M].
# This guess only fails if the number in the closed box is the global minimum or the global maximum
# of all 20 numbers.

# There are 20 numbers in total. The probability that the single closed box contains the
# global minimum is 1/20. The probability that it contains the global maximum is 1/20.
# The total probability of failure is the sum of these two probabilities.
# The general formula for the probability of success when opening k boxes is p = 1 - 2/(k+1).

# Let's calculate this for the optimal k=19.
# We will print each number involved in the final calculation.
one = 1
two = 2
nineteen = 19

# The equation for the probability p is:
# p = 1 - 2 / (19 + 1)
probability = one - two / (nineteen + one)

print("The optimal strategy for Alice is to open 19 boxes.")
print("The maximal probability p of success is given by the equation:")
print(f"p = {one} - {two} / ({nineteen} + {one})")
print(f"p = {one - two / (nineteen + one)}")
print(f"The result is {probability:.1f}, which is equivalent to the fraction 9/10.")
