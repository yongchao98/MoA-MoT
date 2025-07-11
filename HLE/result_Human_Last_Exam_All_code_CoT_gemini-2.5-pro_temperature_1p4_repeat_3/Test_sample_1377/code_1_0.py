import numpy as np
import math

def calculate_optimal_fractions(probabilities, payouts):
    """
    Calculates the Kelly optimal betting fractions for a set of mutually exclusive events.
    
    Args:
        probabilities (list): A list of probabilities for each outcome.
        payouts (list): A list of net payout ratios (e.g., 4 for 4:1) for each outcome.

    Returns:
        np.array: An array containing the optimal fraction of capital to bet on each outcome.
    """
    num_outcomes = len(probabilities)
    p = np.array(probabilities)
    b = np.array(payouts)

    # The optimal fractions f satisfy the system M*f = p/lambda - 1,
    # where M is a matrix derived from the payouts b, and lambda is a
    # constant to be found.
    M = np.diag(b) - np.ones((num_outcomes, num_outcomes)) + np.diag(np.ones(num_outcomes))

    # We can solve for f in terms of lambda: f = M_inv * (p/lambda - 1)
    M_inv = np.linalg.inv(M)
    
    # f = (1/lambda) * (M_inv @ p) - (M_inv @ 1)
    term1 = M_inv @ p
    term2 = M_inv @ np.ones(num_outcomes)

    # To maximize growth, we must minimize lambda subject to constraints:
    # 1. f_i >= 0 for all i
    # 2. sum(f_i) <= 1
    
    # Constraint 2 gives a lower bound for lambda
    sum_term1 = np.sum(term1)
    sum_term2 = np.sum(term2)
    # lambda >= sum_term1 / (1 + sum_term2)
    lower_bound_lambda = sum_term1 / (1 + sum_term2)
    
    # We choose the minimal possible lambda, which is the lower bound.
    optimal_lambda = lower_bound_lambda
    
    # Calculate final fractions with the determined optimal lambda
    fractions = (1 / optimal_lambda) * term1 - term2
    
    # Due to floating point arithmetic, fractions might be slightly off.
    # We can clamp values very close to zero.
    fractions[fractions < 1e-9] = 0
    
    return fractions

def calculate_growth_rate(fractions, probabilities, payouts):
    """
    Calculates the expected log-growth rate for a given betting strategy.
    
    Args:
        fractions (list): The fraction of capital bet on each outcome.
        probabilities (list): The true probabilities of each outcome.
        payouts (list): The net payout ratios for each outcome.
        
    Returns:
        float: The expected log-growth rate.
        list: The list of wealth outcomes for each event.
    """
    total_fraction_bet = np.sum(fractions)
    wealth_outcomes = []
    for i in range(len(fractions)):
        # Wealth if outcome i occurs
        wealth = (1 - total_fraction_bet) + fractions[i] * (payouts[i] + 1)
        wealth_outcomes.append(wealth)
        
    growth_rate = 0
    for i in range(len(fractions)):
        if wealth_outcomes[i] <= 0:
            return -float('inf'), wealth_outcomes # Ruin
        growth_rate += probabilities[i] * math.log(wealth_outcomes[i])
        
    return growth_rate, wealth_outcomes

# --- Main Calculation ---
# 1. Define problem parameters
p_true = [1/2, 1/4, 1/4]
p_believed = [1/4, 1/2, 1/4]
b_payouts = [4, 3, 3]

# 2. Calculate W* (Optimal Growth Rate)
# Bettor uses true probabilities to determine fractions
f_optimal = calculate_optimal_fractions(p_true, b_payouts)
W_star, W_star_outcomes = calculate_growth_rate(f_optimal, p_true, b_payouts)

# 3. Calculate W (Actual Growth Rate)
# Bettor uses believed probabilities to determine fractions
f_mistaken = calculate_optimal_fractions(p_believed, b_payouts)
# The growth rate is evaluated against what actually happens (true probabilities)
W_actual, W_actual_outcomes = calculate_growth_rate(f_mistaken, p_true, b_payouts)

# 4. Calculate the difference
difference = W_star - W_actual

# 5. Print the results clearly, showing the equation numbers
print("Calculation of W* (Optimal Growth Rate):")
print(f"Optimal Fractions f* = ({f_optimal[0]:.5f}, {f_optimal[1]:.5f}, {f_optimal[2]:.5f})")
print(f"W* = {p_true[0]}*log({W_star_outcomes[0]:.5f}) + {p_true[1]}*log({W_star_outcomes[1]:.5f}) + {p_true[2]}*log({W_star_outcomes[2]:.5f})")
print(f"W* = {W_star:.5f}\n")

print("Calculation of W (Actual Growth Rate):")
print(f"Mistaken Fractions f = ({f_mistaken[0]:.5f}, {f_mistaken[1]:.5f}, {f_mistaken[2]:.5f})")
print(f"W = {p_true[0]}*log({W_actual_outcomes[0]:.5f}) + {p_true[1]}*log({W_actual_outcomes[1]:.5f}) + {p_true[2]}*log({W_actual_outcomes[2]:.5f})")
print(f"W = {W_actual:.5f}\n")

print("Final Result:")
print(f"W* - W = {W_star:.5f} - {W_actual:.5f} = {difference:.5f}")

<<<0.22735>>>