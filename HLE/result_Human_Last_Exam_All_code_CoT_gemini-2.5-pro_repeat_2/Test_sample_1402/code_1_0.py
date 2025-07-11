# In a single game:
# The probability of Theo winning is 1/6.
# Therefore, the probability of Theo NOT winning is 1 - 1/6 = 5/6.

p_not_win_numerator = 5
p_not_win_denominator = 6
num_games = 4

# The problem asks for the probability that Theo's first win occurs after at least five games.
# This is equivalent to Theo not winning the first four games.
# We calculate this as (5/6)^4.

final_numerator = p_not_win_numerator ** num_games
final_denominator = p_not_win_denominator ** num_games

print("The probability of Theo not winning a single game is {}/{}.".format(p_not_win_numerator, p_not_win_denominator))
print("The number of games he must not win in a row is {}.".format(num_games))
print("\nThe final probability is calculated by the equation: ({}/{})^{} = {}/{}".format(
    p_not_win_numerator, p_not_win_denominator, num_games, final_numerator, final_denominator
))