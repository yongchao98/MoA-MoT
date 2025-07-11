import locale

# Use locale to format numbers with commas for readability
locale.setlocale(locale.LC_ALL, '')

# --- Define Accredited Investor Thresholds for Ontario ---
INDIVIDUAL_NET_FINANCIAL_ASSETS_THRESHOLD = 1_000_000
INDIVIDUAL_NET_ASSETS_THRESHOLD = 5_000_000
JOINT_INCOME_THRESHOLD = 300_000
ENTITY_NET_ASSETS_THRESHOLD = 5_000_000

def f(n):
    """Helper function to format numbers as currency."""
    return locale.format_string("$%.2f", n, grouping=True)

print("--- Accredited Investor Status Analysis ---\n")

# --- Option A: Limited Partnership ---
print("--- Evaluating Option A: Limited Partnership ---")
liam_net_financial_assets = 5_400_000
jack_total_assets = 18_000_000
jack_total_liabilities = 5_000_000
ace_net_financial_assets = 25_000_000
gift_per_partner = 2_000_000

# Check each partner's status
liam_is_ai = liam_net_financial_assets > INDIVIDUAL_NET_FINANCIAL_ASSETS_THRESHOLD
print(f"Liam's net financial assets of {f(liam_net_financial_assets)} are greater than {f(INDIVIDUAL_NET_FINANCIAL_ASSETS_THRESHOLD)}. Liam is an AI: {liam_is_ai}")

jack_net_assets = jack_total_assets - jack_total_liabilities
jack_is_ai = jack_net_assets > INDIVIDUAL_NET_ASSETS_THRESHOLD
print(f"Jack's net assets ({f(jack_total_assets)} - {f(jack_total_liabilities)} = {f(jack_net_assets)}) are greater than {f(INDIVIDUAL_NET_ASSETS_THRESHOLD)}. Jack is an AI: {jack_is_ai}")

ace_is_ai = ace_net_financial_assets > INDIVIDUAL_NET_FINANCIAL_ASSETS_THRESHOLD
print(f"Ace's net financial assets of {f(ace_net_financial_assets)} are greater than {f(INDIVIDUAL_NET_FINANCIAL_ASSETS_THRESHOLD)}. Ace is an AI: {ace_is_ai}")

gp_corp_assets = gift_per_partner * 3
gp_corp_is_ai = gp_corp_assets > ENTITY_NET_ASSETS_THRESHOLD
print(f"The General Partner corporation's assets ({f(gift_per_partner)} * 3 = {f(gp_corp_assets)}) are greater than {f(ENTITY_NET_ASSETS_THRESHOLD)}. The GP is an AI: {gp_corp_is_ai}")

option_a_is_ai = all([liam_is_ai, jack_is_ai, ace_is_ai, gp_corp_is_ai])
print(f"\nConclusion for A: All owners are accredited investors. Therefore, the limited partnership IS an Accredited Investor. Status: {option_a_is_ai}\n")


# --- Option B: Individual (Joint Income) ---
print("--- Evaluating Option B: Individual (Joint Income) ---")
individual_income_2019 = 150_000
spouse_income_2019 = 170_000
individual_income_2020 = 175_000
spouse_income_2020 = 175_000

joint_income_2019 = individual_income_2019 + spouse_income_2019
joint_income_2020 = individual_income_2020 + spouse_income_2020

income_test_2019 = joint_income_2019 > JOINT_INCOME_THRESHOLD
print(f"2019 joint income: {f(individual_income_2019)} + {f(spouse_income_2019)} = {f(joint_income_2019)}. This is greater than {f(JOINT_INCOME_THRESHOLD)}: {income_test_2019}")

income_test_2020 = joint_income_2020 > JOINT_INCOME_THRESHOLD
print(f"2020 joint income: {f(individual_income_2020)} + {f(spouse_income_2020)} = {f(joint_income_2020)}. This is greater than {f(JOINT_INCOME_THRESHOLD)}: {income_test_2020}")

option_b_is_ai = income_test_2019 and income_test_2020
print(f"\nConclusion for B: Joint income exceeded the threshold for both years. Therefore, the individual IS an Accredited Investor. Status: {option_b_is_ai}\n")


# --- Option C: Individual (Joint Net Assets) ---
print("--- Evaluating Option C: Individual (Joint Net Assets) ---")
total_assets = 6_000_000
total_liabilities = 1_000_000

net_assets = total_assets - total_liabilities
option_c_is_ai = net_assets > INDIVIDUAL_NET_ASSETS_THRESHOLD
print(f"Net Assets calculation: {f(total_assets)} - {f(total_liabilities)} = {f(net_assets)}")
print(f"The rule requires net assets to be GREATER THAN {f(INDIVIDUAL_NET_ASSETS_THRESHOLD)}.")
print(f"Is {f(net_assets)} > {f(INDIVIDUAL_NET_ASSETS_THRESHOLD)}? {option_c_is_ai}")
print(f"\nConclusion for C: Net assets are exactly equal to, but do not exceed, the threshold. Therefore, the individual is NOT an Accredited Investor. Status: {option_c_is_ai}\n")


# --- Option D: Corporation ---
print("--- Evaluating Option D: Corporation ---")
jose_nfa = 100_000_000
transfer_pct = 0.10
corp_assets_from_jose = jose_nfa * transfer_pct

# First, check look-through (does not apply if even one owner is not an AI)
james_income_2019 = 75_000
james_spouse_income_2019 = 210_000
james_joint_income = james_income_2019 + james_spouse_income_2019
james_is_ai = james_joint_income > JOINT_INCOME_THRESHOLD
print(f"Checking James (shareholder): Joint income is {f(james_joint_income)}, which is not > {f(JOINT_INCOME_THRESHOLD)}. James is not an AI: {not james_is_ai}. Look-through test fails.")

# Second, check the corporation's own net assets
option_d_is_ai = corp_assets_from_jose > ENTITY_NET_ASSETS_THRESHOLD
print(f"Checking corporation assets: Jose transferred {f(jose_nfa)} * {transfer_pct:.0%} = {f(corp_assets_from_jose)}")
print(f"The corporation's net assets of {f(corp_assets_from_jose)} are greater than {f(ENTITY_NET_ASSETS_THRESHOLD)}: {option_d_is_ai}")
print(f"\nConclusion for D: Although the look-through test fails, the corporation qualifies on its own. Therefore, the corporation IS an Accredited Investor. Status: {option_d_is_ai}\n")


# --- Option E: Corporation ---
print("--- Evaluating Option E: Corporation ---")
corp_assets = 5_500_000
corp_liabilities = 300_000
corp_net_assets = corp_assets - corp_liabilities

# As established in the prompt, neither Alex nor Bob are accredited investors, so the look-through test fails.
# Check the corporation's own net assets
option_e_is_ai = corp_net_assets > ENTITY_NET_ASSETS_THRESHOLD
print(f"Checking corporation assets: {f(corp_assets)} - {f(corp_liabilities)} = {f(corp_net_assets)}")
print(f"The corporation's net assets of {f(corp_net_assets)} are greater than {f(ENTITY_NET_ASSETS_THRESHOLD)}: {option_e_is_ai}")
print(f"\nConclusion for E: Although the shareholders are not AIs, the corporation qualifies on its own. Therefore, the corporation IS an Accredited Investor. Status: {option_e_is_ai}\n")


print("--- FINAL ANSWER ---")
print("The option that would NOT be classified as an Accredited Investor is C.")
print("<<<C>>>")