import math

class FlightRevenueOptimizer:
    """
    Solves the airline revenue management problem using dynamic programming.
    """
    def __init__(self):
        # Customer demand distributions (number_of_requests, probability)
        self.demand_c1 = [(0, 0.25), (1, 0.5), (2, 0.25)]
        self.demand_c2 = [(0, 0.25), (1, 0.5), (2, 0.25)]
        self.demand_c2_week1 = [(0, 1.0)] # No Class 2 customers in the first week

        # Constants from the problem description
        self.max_seats = 10
        self.total_days = 14
        self.cheap_ticket_price = 100
        self.expensive_ticket_price = 200
        self.c2_expensive_buy_prob = 0.5
        
        # Memoization tables for dynamic programming and combinatorics
        self.memo_dp = {}
        self.memo_nCr = {}

    def nCr(self, n, r):
        """Helper function for combinations nCr, with memoization."""
        if r < 0 or r > n:
            return 0
        state = (n, r)
        if state in self.memo_nCr:
            return self.memo_nCr[state]
        
        # Using math.comb for efficiency and correctness
        result = math.comb(n, r)
        self.memo_nCr[state] = result
        return result

    def calculate_expected_revenue(self, day, cheap_avail, exp_avail):
        """
        Recursive function to calculate expected revenue using dynamic programming.
        Returns a tuple: (total_expected_revenue, cheap_ticket_revenue, expensive_ticket_revenue)
        """
        # Base Case: If the sales period is over, no more revenue can be generated.
        if day > self.total_days:
            return (0, 0, 0)
        
        # Memoization: Return stored result if this state has been computed before.
        state = (day, cheap_avail, exp_avail)
        if state in self.memo_dp:
            return self.memo_dp[state]

        # Determine Class 2 customer demand based on the week.
        current_demand_c2 = self.demand_c2 if day > 7 else self.demand_c2_week1

        # Initialize accumulators for expected revenues for the current state.
        total_rev_acc = 0
        cheap_rev_acc = 0
        exp_rev_acc = 0
        
        # Iterate over all possible demand scenarios for Class 1 and Class 2 customers.
        for n1, p1 in self.demand_c1:
            for n2, p2 in current_demand_c2:
                prob_scenario = p1 * p2
                
                # --- Step 1: Handle cheap ticket sales ---
                # Class 2 has priority for cheap tickets.
                c_sold_to_2 = min(n2, cheap_avail)
                c_rem_for_c1 = cheap_avail - c_sold_to_2
                
                # Class 1 buys the remaining cheap tickets.
                c_sold_to_1 = min(n1, c_rem_for_c1)
                
                rev_cheap_today = (c_sold_to_1 + c_sold_to_2) * self.cheap_ticket_price
                c_final_avail = cheap_avail - c_sold_to_1 - c_sold_to_2

                # --- Step 2: Handle expensive ticket sales ---
                n2_exp_demand = n2 - c_sold_to_2
                
                # Accumulators for the expectation over the binomial distribution of C2 choices.
                inner_total_exp, inner_cheap_exp, inner_exp_exp = 0, 0, 0

                # If C2 customers might buy expensive tickets and such tickets are available.
                if n2_exp_demand > 0 and exp_avail > 0:
                    N = n2_exp_demand
                    # Iterate over all possible outcomes (k) of N customers deciding to buy.
                    for k in range(N + 1):
                        # Probability that exactly k customers decide to buy (from Binomial distribution).
                        prob_k = self.nCr(N, k) * (self.c2_expensive_buy_prob ** N)

                        e_sold = min(k, exp_avail)
                        rev_exp_today_k = e_sold * self.expensive_ticket_price
                        exp_final_avail = exp_avail - e_sold
                        
                        future_revs = self.calculate_expected_revenue(day + 1, c_final_avail, exp_final_avail)
                        
                        inner_total_exp += prob_k * (rev_exp_today_k + future_revs[0])
                        inner_cheap_exp += prob_k * future_revs[1]
                        inner_exp_exp += prob_k * (rev_exp_today_k + future_revs[2])
                else:
                    future_revs = self.calculate_expected_revenue(day + 1, c_final_avail, exp_avail)
                    inner_total_exp, inner_cheap_exp, inner_exp_exp = future_revs

                # --- Step 3: Aggregate results for the current state ---
                total_rev_acc += prob_scenario * (rev_cheap_today + inner_total_exp)
                cheap_rev_acc += prob_scenario * (rev_cheap_today + inner_cheap_exp)
                exp_rev_acc += prob_scenario * inner_exp_exp

        # Memoize and return the calculated expected revenues for the state.
        self.memo_dp[state] = (total_rev_acc, cheap_rev_acc, exp_rev_acc)
        return self.memo_dp[state]

    def find_max_revenue(self):
        """
        Iterates through all possible booking policies to find the one that
        maximizes total expected revenue.
        """
        best_policy = {"b": -1, "total_rev": 0, "cheap_rev": 0, "exp_rev": 0}

        # A policy 'b' is the number of cheap tickets offered. Test all possibilities.
        for b in range(self.max_seats + 1):
            self.memo_dp.clear() # Reset memoization for each new policy.
            exp_revs = self.calculate_expected_revenue(1, b, self.max_seats - b)
            
            if exp_revs[0] > best_policy["total_rev"]:
                best_policy = {
                    "b": b,
                    "total_rev": exp_revs[0],
                    "cheap_rev": exp_revs[1],
                    "exp_rev": exp_revs[2]
                }
        
        print(f"The optimal policy is to offer {best_policy['b']} cheap tickets and {self.max_seats - best_policy['b']} expensive tickets.")
        
        cheap_part = best_policy['cheap_rev']
        exp_part = best_policy['exp_rev']
        total_part = best_policy['total_rev']

        print(f"Expected Revenue from Cheap Tickets = {cheap_part:.4f}")
        print(f"Expected Revenue from Expensive Tickets = {exp_part:.4f}")
        print(f"Maximum Expected Total Revenue = {cheap_part:.4f} + {exp_part:.4f} = {total_part:.4f}")

        return total_part


# Execute the solver
solver = FlightRevenueOptimizer()
max_revenue = solver.find_max_revenue()
print(f"<<<{max_revenue:.4f}>>>")