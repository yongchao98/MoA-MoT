import random
import numpy as np

# --- Configuration ---
LOAN_COUNT = 1000
LOAN_PRINCIPAL = 300000.0

class Mortgage:
    """Represents a single mortgage with its key risk attributes."""
    def __init__(self, principal, fico, ltv):
        self.principal = principal
        self.fico = fico
        self.ltv = ltv
        self.is_defaulted = False

class RMBS:
    """Represents a Residential Mortgage-Backed Security holding a pool of mortgages."""
    def __init__(self, name, mortgages):
        self.name = name
        self.mortgages = mortgages

    def get_value(self):
        """Calculates the current value based on non-defaulted loans."""
        return sum(m.principal for m in self.mortgages if not m.is_defaulted)

    def get_avg_fico(self):
        """Calculates the average FICO score of the pool."""
        return sum(m.fico for m in self.mortgages) / len(self.mortgages)

    def get_avg_ltv(self):
        """Calculates the average Loan-to-Value ratio of the pool."""
        return sum(m.ltv for m in self.mortgages) / len(self.mortgages)

    def get_default_rate(self):
        """Calculates the percentage of defaulted loans."""
        defaulted_count = sum(1 for m in self.mortgages if m.is_defaulted)
        return defaulted_count / len(self.mortgages)


def create_mortgage_pool(originator_quality, count, principal):
    """Creates a pool of mortgages based on the originator's quality."""
    mortgages = []
    for _ in range(count):
        if originator_quality == "High":
            # High-quality originators have strict standards
            fico = int(random.normalvariate(760, 25))
            ltv = max(0, min(1, random.normalvariate(0.70, 0.05)))
        else: # Low quality
            # Low-quality originators have lax standards ("subprime")
            fico = int(random.normalvariate(620, 40))
            ltv = max(0, min(1, random.normalvariate(0.95, 0.07)))
        mortgages.append(Mortgage(principal, fico, ltv))
    return mortgages

def run_market_stress_test(mortgage_pool):
    """Simulates a market downturn, causing defaults based on FICO and LTV."""
    for mortgage in mortgage_pool:
        # Simplified default probability model
        # Risk increases significantly for FICO < 660 and LTV > 0.90
        risk_score = 0.0
        if mortgage.fico < 660:
            risk_score += 0.25 # Base risk for subprime FICO
        elif mortgage.fico < 720:
            risk_score += 0.05 # Base risk for average FICO

        if mortgage.ltv > 0.90:
            risk_score += 0.30 # Additional risk for high LTV
        elif mortgage.ltv > 0.80:
            risk_score += 0.10 # Additional risk for moderate LTV

        if random.random() < risk_score:
            mortgage.is_defaulted = True

# --- Main Simulation ---

print("--- Step 1: Origination ---")
print("Two RMBS pools are created, each with 1000 mortgages of $300,000.")
print("The only difference is the quality of the loan originator.")

# Create pools from two different types of originators
high_quality_pool = create_mortgage_pool("High", LOAN_COUNT, LOAN_PRINCIPAL)
low_quality_pool = create_mortgage_pool("Low", LOAN_COUNT, LOAN_PRINCIPAL)

# Package them into RMBS
rmbs_a = RMBS("RMBS from High-Quality Originator", high_quality_pool)
rmbs_b = RMBS("RMBS from Low-Quality Originator", low_quality_pool)

initial_value = LOAN_COUNT * LOAN_PRINCIPAL

print("\n--- Step 2: Initial Pool Characteristics ---")
print(f"'{rmbs_a.name}':")
print(f"  Initial Value: ${initial_value:,.0f}")
print(f"  Average FICO Score: {rmbs_a.get_avg_fico():.0f}")
print(f"  Average LTV Ratio: {rmbs_a.get_avg_ltv():.2%}")

print(f"\n'{rmbs_b.name}':")
print(f"  Initial Value: ${initial_value:,.0f}")
print(f"  Average FICO Score: {rmbs_b.get_avg_fico():.0f}")
print(f"  Average LTV Ratio: {rmbs_b.get_avg_ltv():.2%}")
print("\nNote how the Low-Quality originator created loans with lower FICO and higher LTV.")


print("\n--- Step 3: Market Stress Test ---")
print("A market downturn is simulated (e.g., housing prices fall, interest rates rise).")
print("Defaults are calculated based on the underlying loan quality (FICO, LTV).")

run_market_stress_test(rmbs_a.mortgages)
run_market_stress_test(rmbs_b.mortgages)

print("\n--- Step 4: Post-Stress Results ---")
print("The simulation shows the devastating impact on the RMBS from the low-quality originator.")

value_a = rmbs_a.get_value()
loss_a = initial_value - value_a
loss_pct_a = (loss_a / initial_value)
print(f"\n'{rmbs_a.name}':")
print(f"  Post-Stress Value: ${value_a:,.0f}")
print(f"  Default Rate: {rmbs_a.get_default_rate():.2%}")
print(f"  Value Loss: ${loss_a:,.0f} ({loss_pct_a:.2%})")

value_b = rmbs_b.get_value()
loss_b = initial_value - value_b
loss_pct_b = (loss_b / initial_value)
print(f"\n'{rmbs_b.name}':")
print(f"  Post-Stress Value: ${value_b:,.0f}")
print(f"  Default Rate: {rmbs_b.get_default_rate():.2%}")
print(f"  Value Loss: ${loss_b:,.0f} ({loss_pct_b:.2%})")

print("\n--- Conclusion ---")
print("The simulation demonstrates that the quality of the loan originator is the root cause factor.")
print("Lax underwriting standards directly created a fragile RMBS that collapsed under stress, while the RMBS with properly underwritten loans remained much more stable.")