import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Create a string buffer
output_buffer = io.StringIO()
# Redirect stdout
sys.stdout = output_buffer


def format_currency(value):
    """Formats a number as a currency string."""
    return f"${value:,.2f}"

# Define Accredited Investor thresholds
FINANCIAL_ASSETS_THRESHOLD = 1_000_000
NET_ASSETS_THRESHOLD = 5_000_000
INDIVIDUAL_INCOME_THRESHOLD = 200_000
JOINT_INCOME_THRESHOLD = 300_000
CORPORATION_NET_ASSETS_THRESHOLD = 5_000_000

print("Analyzing each option based on Accredited Investor definitions in Ontario...\n")

# --- Option A: Limited Partnership ---
print("--- Analysis of Option A ---")
liam_net_financial_assets = 5_400_000
jack_total_assets = 18_000_000
jack_total_liabilities = 5_000_000
ace_net_financial_assets = 25_000_000
gp_gift_per_partner = 2_000_000
gp_partner_count = 3

jack_net_assets = jack_total_assets - jack_total_liabilities
gp_corp_net_assets = gp_gift_per_partner * gp_partner_count

liam_is_ai = liam_net_financial_assets >= FINANCIAL_ASSETS_THRESHOLD
jack_is_ai = jack_net_assets >= NET_ASSETS_THRESHOLD
ace_is_ai = ace_net_financial_assets >= FINANCIAL_ASSETS_THRESHOLD
gp_corp_is_ai = gp_corp_net_assets >= CORPORATION_NET_ASSETS_THRESHOLD

partnership_is_ai_via_lookthrough = liam_is_ai and jack_is_ai and ace_is_ai and gp_corp_is_ai

print("A limited partnership can qualify if all its partners are accredited investors.")
print(f"Liam's net financial assets are {format_currency(liam_net_financial_assets)}, which is >= {format_currency(FINANCIAL_ASSETS_THRESHOLD)}. Liam is an AI.")
print(f"Jack's net assets are {format_currency(jack_total_assets)} - {format_currency(jack_total_liabilities)} = {format_currency(jack_net_assets)}, which is >= {format_currency(NET_ASSETS_THRESHOLD)}. Jack is an AI.")
print(f"Ace's net financial assets are {format_currency(ace_net_financial_assets)}, which is >= {format_currency(FINANCIAL_ASSETS_THRESHOLD)}. Ace is an AI.")
print(f"The General Partner corporation's net assets are {gp_partner_count} * {format_currency(gp_gift_per_partner)} = {format_currency(gp_corp_net_assets)}, which is >= {format_currency(CORPORATION_NET_ASSETS_THRESHOLD)}. The GP is an AI.")
print("Conclusion: Since all partners are accredited investors, the limited partnership in Option A is an Accredited Investor.\n")


# --- Option B: Individual (Joint Income) ---
print("--- Analysis of Option B ---")
individual_income_2019 = 150_000
spouse_income_2019 = 170_000
individual_income_2020 = 175_000
spouse_income_2020 = 175_000

joint_income_2019 = individual_income_2019 + spouse_income_2019
joint_income_2020 = individual_income_2020 + spouse_income_2020

is_ai_2019 = joint_income_2019 > JOINT_INCOME_THRESHOLD
is_ai_2020 = joint_income_2020 > JOINT_INCOME_THRESHOLD

individual_is_ai = is_ai_2019 and is_ai_2020

print("An individual can qualify if their joint income with a spouse exceeded $300,000 in the last two years.")
print(f"2019 Joint Income: {format_currency(individual_income_2019)} + {format_currency(spouse_income_2019)} = {format_currency(joint_income_2019)}, which is > {format_currency(JOINT_INCOME_THRESHOLD)}.")
print(f"2020 Joint Income: {format_currency(individual_income_2020)} + {format_currency(spouse_income_2020)} = {format_currency(joint_income_2020)}, which is > {format_currency(JOINT_INCOME_THRESHOLD)}.")
print("Conclusion: Since the income test is met for both years, the individual in Option B is an Accredited Investor.\n")


# --- Option C: Individual (Joint Net Assets) ---
print("--- Analysis of Option C ---")
couple_total_assets = 6_000_000
couple_total_liabilities = 1_000_000

couple_net_assets = couple_total_assets - couple_total_liabilities
individual_is_ai_c = couple_net_assets >= NET_ASSETS_THRESHOLD

print("An individual can qualify if their joint net assets with a spouse are at least $5,000,000.")
print(f"Joint Net Assets: {format_currency(couple_total_assets)} - {format_currency(couple_total_liabilities)} = {format_currency(couple_net_assets)}.")
print(f"This value of {format_currency(couple_net_assets)} meets the threshold of {format_currency(NET_ASSETS_THRESHOLD)}.")
print("Conclusion: The individual in Option C is an Accredited Investor.\n")


# --- Option D: Corporation ---
print("--- Analysis of Option D ---")
jose_net_financial_assets = 100_000_000
transfer_percentage = 0.10
james_income = 75_000
james_spouse_income = 210_000

corp_assets_from_jose = jose_net_financial_assets * transfer_percentage
corp_d_is_ai_on_assets = corp_assets_from_jose >= CORPORATION_NET_ASSETS_THRESHOLD

james_joint_income = james_income + james_spouse_income
james_is_ai = james_joint_income > JOINT_INCOME_THRESHOLD

print("A corporation can qualify if its net assets are at least $5,000,000.")
print(f"Jose transferred {transfer_percentage:.0%} of his assets to the corporation: {format_currency(jose_net_financial_assets)} * {transfer_percentage:.0%} = {format_currency(corp_assets_from_jose)}.")
print(f"The corporation's net assets of {format_currency(corp_assets_from_jose)} are greater than the {format_currency(CORPORATION_NET_ASSETS_THRESHOLD)} threshold.")
print("Even though one shareholder, James, is not an AI (joint income of {format_currency(james_income)} + {format_currency(james_spouse_income)} = {format_currency(james_joint_income)} is less than {format_currency(JOINT_INCOME_THRESHOLD)}), the corporation qualifies on its own assets.")
print("Conclusion: The corporation in Option D is an Accredited Investor.\n")


# --- Option E: Corporation ---
print("--- Analysis of Option E ---")
corp_e_gic = 5_500_000
corp_e_liabilities = 300_000
alex_net_financial_assets = 900_000
alex_net_assets = 3_000_000
bob_income_2020 = 41_000
bob_financial_assets = 75_000

corp_e_net_assets = corp_e_gic - corp_e_liabilities
corp_e_is_ai_on_assets = (corp_e_gic - corp_e_liabilities) >= CORPORATION_NET_ASSETS_THRESHOLD
alex_is_ai = (alex_net_financial_assets >= FINANCIAL_ASSETS_THRESHOLD) or (alex_net_assets >= NET_ASSETS_THRESHOLD)
bob_is_ai = (bob_income_2020 > INDIVIDUAL_INCOME_THRESHOLD) or (bob_financial_assets >= FINANCIAL_ASSETS_THRESHOLD)

print("A corporation can qualify on its own net assets or if all its shareholders are AIs.")
print(f"1. Corporation's own assets: Its minimum net assets are at least {format_currency(corp_e_gic)} - {format_currency(corp_e_liabilities)} = {format_currency(corp_e_net_assets)} (plus the value of properties).")
print(f"This amount ({format_currency(corp_e_net_assets)}) is greater than the {format_currency(CORPORATION_NET_ASSETS_THRESHOLD)} threshold. On this basis alone, it seems to qualify.")
print("\n2. Shareholders' status ('Look-through' test):")
print(f"   - Alex's net financial assets ({format_currency(alex_net_financial_assets)}) and net assets ({format_currency(alex_net_assets)}) are both below the required thresholds. Alex is NOT an AI.")
print(f"   - Bob's income and assets are far below the thresholds. Bob is NOT an AI.")
print("\nConclusion: Although the corporation's assets appear to meet the test, it is owned entirely by individuals who are NOT accredited investors. Securities regulations include anti-avoidance provisions to prevent non-accredited individuals from pooling assets simply to overcome the thresholds. Because this entity is a holding company for two non-accredited investors, it would likely NOT be classified as an Accredited Investor by regulators, as it defeats the purpose of the rule.")
print("Therefore, Option E represents the entity that would not be classified as an Accredited Investor.\n")

final_answer = "E"

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output = output_buffer.getvalue()

# Print the captured output
print(output)
print(f"<<<{final_answer}>>>")