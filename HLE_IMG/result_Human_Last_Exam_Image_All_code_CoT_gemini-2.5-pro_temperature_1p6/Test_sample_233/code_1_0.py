# The task is to identify the Ingvar runestone from the provided image fragment.
# Based on the analysis of the runes, key text fragments were identified:
# 1. 'þuriþ' (the name Þuríðr/Thorid)
# 2. 'stain' (stæin/stone)
# 3. 'aftir' (æftiR/in memory of)
#
# A search through the corpus of Ingvar runestones for these fragments
# leads to the identification of the Strängnäs Cathedral runestone fragment.

# The ID of a Swedish runestone consists of a provincial code and a number.
province_code = "Sö"  # Sö stands for the province of Södermanland.
runestone_number = 279

# The script will now print the components and the final ID.
print(f"The runestone is located in the province of Södermanland, indicated by the code: '{province_code}'")
print(f"The number assigned to this specific runestone in the catalog is: {runestone_number}")
print("-" * 20)
print(f"The full ID of the Ingvar runestone is: {province_code} {runestone_number}")
