# Constants for Accredited Investor (AI) tests in Canadian dollars
FINANCIAL_ASSETS_TEST = 1000000.00
NET_ASSET_TEST = 5000000.00
INCOME_TEST = 200000.00
JOINT_INCOME_TEST = 300000.00
COMPANY_NET_ASSET_TEST = 5000000.00

print("Analyzing each option based on the Accredited Investor definitions in Ontario.")

# --- Option A: The Limited Partnership ---
print("\n--- Analysis of Option A ---")
liam_net_financial_assets = 5400000.00
jack_total_assets = 18000000.00
jack_total_liabilities = 5000000.00
ace_net_financial_assets = 25000000.00
gp_gift_total = 2000000.00 * 3

# Check Limited Partners
liam_is_ai = liam_net_financial_assets > FINANCIAL_ASSETS_TEST
jack_net_assets = jack_total_assets - jack_total_liabilities
jack_is_ai = jack_net_assets >= NET_ASSET_TEST
ace_is_ai = ace_net_financial_assets > FINANCIAL_ASSETS_TEST
gp_is_ai = gp_gift_total >= COMPANY_NET_ASSET_TEST

all_owners_are_ai = liam_is_ai and jack_is_ai and ace_is_ai and gp_is_ai

print(f"Liam's net financial assets: ${liam_net_financial_assets:,.2f}. Meets >${FINANCIAL_ASSETS_TEST:,.2f} test. Liam is an AI.")
print(f"Jack's net assets: ${jack_total_assets:,.2f} - ${jack_total_liabilities:,.2f} = ${jack_net_assets:,.2f}. Meets >=${NET_ASSET_TEST:,.2f} test. Jack is an AI.")
print(f"Ace's net financial assets: ${ace_net_financial_assets:,.2f}. Meets >${FINANCIAL_ASSETS_TEST:,.2f} test. Ace is an AI.")
print(f"General Partner's net assets: ${gp_gift_total:,.2f}. Meets >=${COMPANY_NET_ASSET_TEST:,.2f} test. The GP is an AI.")
print(f"Conclusion for A: The Limited Partnership qualifies because all of its owners are Accredited Investors. It IS an AI.")


# --- Option B: The Individual (Joint Income) ---
print("\n--- Analysis of Option B ---")
individual_income_2019 = 150000.00
individual_income_2020 = 175000.00
spouse_income_2019 = 170000.00
spouse_income_2020 = 175000.00

joint_income_2019 = individual_income_2019 + spouse_income_2019
joint_income_2020 = individual_income_2020 + spouse_income_2020

print(f"Combined income in 2019: ${individual_income_2019:,.2f} + ${spouse_income_2019:,.2f} = ${joint_income_2019:,.2f}")
print(f"Combined income in 2020: ${individual_income_2020:,.2f} + ${spouse_income_2020:,.2f} = ${joint_income_2020:,.2f}")
print(f"Conclusion for B: The individual qualifies because joint income exceeded ${JOINT_INCOME_TEST:,.2f} in the last two years. It IS an AI.")


# --- Option C: The Individual (Joint Net Assets) ---
print("\n--- Analysis of Option C ---")
total_assets = 6000000.00
total_liabilities = 1000000.00

net_assets = total_assets - total_liabilities

print(f"Combined net assets: ${total_assets:,.2f} - ${total_liabilities:,.2f} = ${net_assets:,.2f}")
print(f"Conclusion for C: The individual qualifies because joint net assets are at least ${NET_ASSET_TEST:,.2f}. It IS an AI.")


# --- Option D: The Corporation (Owners & Assets) ---
print("\n--- Analysis of Option D ---")
jose_net_financial_assets = 100000000.00
james_income_2019 = 75000.00
james_spouse_income = 210000.00
corp_asset_transfer = jose_net_financial_assets * 0.10

james_joint_income = james_income_2019 + james_spouse_income
james_is_ai = james_joint_income > JOINT_INCOME_TEST

print(f"Jose has net financial assets of ${jose_net_financial_assets:,.2f} and is an AI.")
print(f"James's joint income is ${james_income_2019:,.2f} + ${james_spouse_income:,.2f} = ${james_joint_income:,.2f}, which is less than ${JOINT_INCOME_TEST:,.2f}. James is not an AI.")
print(f"The corporation's net assets are ${corp_asset_transfer:,.2f}, which is >= ${COMPANY_NET_ASSET_TEST:,.2f}. This qualifies it as an AI.")
print("Also, a majority of shares (99%) are owned by Jose, an AI, which is another way it qualifies.")
print("Conclusion for D: The corporation qualifies as an AI. It IS an AI.")


# --- Option E: The Corporation (Mixed Owners and Assets) ---
print("\n--- Analysis of Option E ---")
corp_assets_gic = 5500000.00
corp_liabilities = 300000.00
alex_net_financial_assets = 900000.00
alex_net_assets = 3000000.00

alex_is_ai = (alex_net_financial_assets > FINANCIAL_ASSETS_TEST) or (alex_net_assets >= NET_ASSET_TEST)
# Bob's financials clearly do not meet any test.
bob_is_ai = False

print(f"Alex's net financial assets are ${alex_net_financial_assets:,.2f} (<${FINANCIAL_ASSETS_TEST:,.2f}) and net assets are ${alex_net_assets:,.2f} (<${NET_ASSET_TEST:,.2f}). Alex is not an AI.")
print("Bob's income and assets are too low. Bob is not an AI.")
print("Since neither owner is an AI, the corporation cannot qualify based on its ownership.")
corp_net_assets = corp_assets_gic - corp_liabilities
print(f"Checking the corporation's net assets: ${corp_assets_gic:,.2f} - ${corp_liabilities:,.2f} = ${corp_net_assets:,.2f}")
print(f"This value (${corp_net_assets:,.2f}) is greater than the ${COMPANY_NET_ASSET_TEST:,.2f} threshold.")
print("However, an anti-avoidance rule prevents a company from being used by non-accredited investors to pool assets solely to meet the $5,000,000 net asset test.")
print("Conclusion for E: Because its owners are not AIs and it appears to be a vehicle for pooling assets, the corporation would NOT be classified as an Accredited Investor.")