import numpy as np

def calculate_critical_correlation():
    """
    This function calculates the critical amount of correlation required to balance potentiation and depression.
    
    The steps are:
    1. The evolution of selectivity (difference in weights) is driven by the term (sigma^2 - c), where
       sigma^2 is the variance of the input signals and c is their covariance.
    2. The critical point for balance is when this term is zero, meaning the critical covariance c_crit = sigma^2.
    3. For the specified Poisson process, the variance of the activity is equal to its mean rate (sigma^2 = mu).
    4. The mean rate mu is given by the inverse of the inter-event interval.
    5. Therefore, the critical correlation (covariance c) is equal to the mean rate mu.
    """
    # The inter-event interval for the Poisson process is 150s.
    inter_event_interval = 150

    # The average rate of activation, mu, is the inverse of the interval.
    mu = 1 / inter_event_interval

    # For a Poisson process, the variance (sigma^2) of the activity is equal to its mean rate (mu).
    sigma_squared = mu
    
    # The critical amount of correlation 'c' is the covariance value that balances potentiation and depression.
    # Our analysis shows this occurs when c = sigma^2.
    c_critical = sigma_squared

    # Print the final equation with the calculated numbers
    print(f"The critical amount of correlation (covariance, c) is determined by the condition: c = \u03C3\u00B2")
    print(f"For a Poisson process, \u03C3\u00B2 = \u03BC (variance equals mean rate).")
    print(f"Given an inter-event interval of {inter_event_interval}s, the mean rate \u03BC = 1 / {inter_event_interval}.")
    print(f"Therefore, the critical correlation is c = {c_critical:.4f}")
    
    # The final numerical answer in the required format.
    print(f"\nFinal calculated value: {c_critical}")

calculate_critical_correlation()

# Final answer to be extracted
critical_correlation_value = 1 / 150
print(f'<<<{critical_correlation_value}>>>')